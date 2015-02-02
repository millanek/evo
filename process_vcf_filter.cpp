//
//  process_vcf_filter.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_filter.h"

// Filtering constants
static const int MIN_OVERALL_VARIANT_PHRED_QUAL=30;   // Corresponds to 0.001 probability of error

#define SUBPROGRAM "filter"

static const char *FILTER_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE\n"
"Custom file for VCF filtering, output to STD_OUT\n"
"\n"
"       --help                           display this help and exit\n"
"\n"
"       The filtereing parameters that can be changed are:\n"
"       -d, --overall-max-depth (DEFAULT +Inf)  Maximum read depth allowed at the putative variant site - filters out variants due to\n"
"                                               collapsed repeats in the reference (reads from multiple sites are all mapped to one)\n"
"       -m, --overall-min-depth (DEFAULT 0)     Minimum read depth allowed at the putative variant site - filters out strange regions where very few reads align\n"
"       -c, --min-copies=MIN (DEFAULT 1)        The variant needs to be present in at least MIN copies\n"
"                                               (i.e. setting --min_copies==1 gives all segregating sites including singletons)\n"
"       -s, --min-depth-per-sample=MIN (DEFAULT 0)        Minimum read depth at the variant position for any sample\n"
"                                               (i.e. all samples need to have at least MIN reads at the variant position)\n"
"       --allow-missing                         Allow missing genotypes (./.) where read depth < min-depth-per-sample\n"
"       --keep-triallelic                       Do not filter out sites that are not biallelic\n"
"       --max-het-individuals=N (DEFAULT +Inf)  Filter out sites where more than N individuals are heterozygous\n"
"\n"
"       The filtereing parameters that are hard coded are:\n"
"       MIN_OVERALL_VARIANT_PHRED_QUAL=30   The minimum accepted phred scaled [10*log_10*p(variant)] variant quality\n"
"\n"
"       You can also just output statistics on the first 2,000,000 variants:\n"
"       --stats                                 Get stats to help with selecting parameters for filtering\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_HELP = 1, OPT_TRIALLELIC, OPT_ALLOW_MISSING_GENOTYPES, OPT_MAX_NUM_HET, OPT_STATS };

static const char* shortopts = "d:c:s:m:";

static const struct option longopts[] = {
    { "overall-max-depth",       required_argument, NULL, 'd' },
    { "overall-min-depth",       required_argument, NULL, 'm' },
    { "min-copies",       required_argument, NULL, 'c' },
    { "min-depth-per-sample",       required_argument, NULL, 's' },
    { "help",   no_argument, NULL, OPT_HELP },
    { "max-het-individuals", required_argument, NULL, OPT_MAX_NUM_HET },
    { "keep-triallelic",   no_argument, NULL, OPT_TRIALLELIC },
    { "allow-missing",   no_argument, NULL, OPT_ALLOW_MISSING_GENOTYPES },
    { "stats",   no_argument, NULL, OPT_STATS },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int min_copies=1;
    static int max_overall_depth = std::numeric_limits<int>::max();
    static int min_overall_depth = 0;
    static int max_het_indiv = std::numeric_limits<int>::max();
    static int min_depth_in_any_individual = 0;
    static bool bBiallelicFilter = true;
    static bool bAllowMissingGenotpyes = false;
    static bool bStats = false;
    static string vcfFile;
}



int filterMain(int argc, char** argv) {
    parseFilterOptions(argc, argv);
    //std::ios_base::openmode mode = std::ios_base::in;
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    std::vector<int> depthsHetFailed;
    std::vector<int> depthsHetPassed;
    std::vector<int> numVariantsPerHetCount;
    std::vector<std::vector<int> > num_indiv_het_vs_depth;
    
    std::ofstream* statsFFile;
    std::ofstream* statsVarDepthFile;
    std::ofstream* statsStrandBiasFile;
    std::ofstream* statsVariantQualityFile;
    
    
    if (opt::bStats) {
        statsFFile = new std::ofstream(fileRoot + ".inbreeding");
        statsVarDepthFile = new std::ofstream(fileRoot + ".varDepth");
        statsStrandBiasFile = new std::ofstream(fileRoot + ".strandBias");
        statsVariantQualityFile = new std::ofstream(fileRoot + ".varQual");
        std::cerr << "Getting statistics for variants from: " << fileName << " to help with setting filtering criteria" << std::endl;
    } else {
        std::cerr << "Filtering variants from: " << fileName << std::endl;
        std::cerr << "Maximum read depth (overall) set to: " << opt::max_overall_depth << std::endl;
        std::cerr << "Minimum read depth (overall) set to: " << opt::min_overall_depth << std::endl;
        std::cerr << "Minimum copies for a variant (e.g. 1 for allowing singletons): " << opt::min_copies << std::endl;
        std::cerr << "Minimum read depth at the variant position: " << opt::min_depth_in_any_individual << std::endl;
        if (opt::bAllowMissingGenotpyes)
            std::cerr << "Genotypes for individuals with lower read depth at the variant position will be set to missing (./.)" << std::endl;
        if (opt::max_het_indiv < std::numeric_limits<int>::max())
            std::cerr << "Filter out sites where more than " << opt::max_het_indiv << " individuals are heterozygous"  << std::endl;
    }
    // Open connection to read from the vcf file
    // std::ifstream* inFile = new std::ifstream(fileName.c_str(), mode);
    std::istream* inFile = createReader(fileName);
    
    
    bool gotChromosomeNumber = false;
    int numChromosomes;
    string line;
    int totalVariantNumber = 0;
    while (getline(*inFile, line)) {
        if (line[0] == '#') 
            std::cout << line << std::endl;  
        else {
            totalVariantNumber++;
            FilterResult result;
            std::vector<std::string> fields = split(line, '\t');
            if (!gotChromosomeNumber) {
                const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
                numChromosomes = (int)numSamples * 2;
                std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
                numVariantsPerHetCount.assign(numSamples + 1, 0);
                gotChromosomeNumber = true;
            }
            result.overallQuality = atoi(fields[5].c_str());
            
            result.counts = getThisVariantCounts(fields);
            
            if (opt::bStats) {
                *statsFFile << result.counts.inbreedingCoefficient << std::endl;
                *statsVarDepthFile << result.counts.overallDepth << std::endl;
                *statsStrandBiasFile << result.counts.FSpval << std::endl;
                *statsVariantQualityFile << result.overallQuality << std::endl;
                if (totalVariantNumber >= 1000000)
                    exit(EXIT_SUCCESS);
                else
                    continue;
            }
            
            if (result.overallQuality < MIN_OVERALL_VARIANT_PHRED_QUAL)
                continue;
        
            
            if (opt::bBiallelicFilter) {
                result.biallelicPassed = testBiallelic(fields[4]);
                if (!result.biallelicPassed) continue;
            } else
                result.biallelicPassed = true;
            
            
            if (result.counts.overallDepth >= opt::min_overall_depth && result.counts.overallDepth <= opt::max_overall_depth) {
                result.overallDepthPassed = true;
            } else
                continue;
            
            // filter out sites where more than MAX_NUM_HET individuals are heterozygous
            int num_hets = 0; bool mnhPassed;
            if (opt::max_het_indiv < std::numeric_limits<int>::max()) {
                for (std::vector<std::vector<int> >::size_type i = 0; i < result.counts.individualsWithVariant.size(); i++) {
                    if (result.counts.individualsWithVariant[i] == 1)
                        num_hets++;
                }
                mnhPassed = (num_hets < opt::max_het_indiv) ? true : false;
            } else {
                mnhPassed = true;
            }
            
            
            if (mnhPassed) {
                if (result.biallelicPassed && result.counts.overall <= (numChromosomes - opt::min_copies) && result.counts.overall >= opt::min_copies) {
    // bool mnhPassed = testMaxNumHet(result, depthsHetFailed, depthsHetPassed, numVariantsPerHetCount, opt::max_het_indiv,num_indiv_het_vs_depth);
                    if (result.counts.minimumDepthInAnIndividual >= opt::min_depth_in_any_individual) {
                        std::cout << line << std::endl;
                    } else if (opt::bAllowMissingGenotpyes) {
                        // Set genotypes to missing for individuals with too low coverage
                        for (int i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                            if (result.counts.depthPerIndividual[i-NUM_NON_GENOTYPE_COLUMNS] < opt::min_depth_in_any_individual) {
                                fields[i][0] = '.'; fields[i][2] = '.';
                            }
                        }
                        print_vector_stream(fields, std::cout,'\t');
                    }
                }
            }
        }
    }
    return 0;
}




void parseFilterOptions(int argc, char** argv) {
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'd': arg >> opt::max_overall_depth; break;
            case 'm': arg >> opt::min_overall_depth; break;
            case 'c': arg >> opt::min_copies; break;
            case 's': arg >> opt::min_depth_in_any_individual; break;
            case '?': die = true; break;
            case OPT_TRIALLELIC: opt::bBiallelicFilter = false; break;
            case OPT_ALLOW_MISSING_GENOTYPES: opt::bAllowMissingGenotpyes = true; break;
            case OPT_MAX_NUM_HET: arg >> opt::max_het_indiv; break;
            case OPT_STATS: opt::bStats = true; break;
            case OPT_HELP:
                std::cout << FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 1) 
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << FILTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    
}
