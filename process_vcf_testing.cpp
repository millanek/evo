//
//  process_vcf_testing.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 05/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_testing.h"

// Filtering constants
static const int MIN_OVERALL_VARIANT_PHRED_QUAL=30;   // Corresponds to 0.001 probability of error

#define SUBPROGRAM "test"

static const char *TESTING_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE\n"
"Custom file for VCF filtering, output to STD_OUT\n"
"\n"
"       --help                           display this help and exit\n"
"\n"
"       The filtereing parameters that can be changed are:\n"
"       -d, --overall-max-depth (DEFAULT +Inf)  Maximum read depth allowed at the putative variant site - fisters out variants due to\n"
"                                               collapsed repeats in the reference (reads from multiple sites are all mapped to one)\n"
"       -c, --min-copies=MIN (DEFAULT 1)        The variant needs to be present in at least MIN copies\n"
"                                               (i.e. setting --min_copies==1 gives all segregating sites including singletons)\n"
"       -s, --min-depth-per-sample=MIN (DEFAULT 3)        Minimum read depth at the variant position for any sample\n"
"                                               (i.e. all samples neet to have at least MIN reads at the variant position)\n"
"       --max-het-individuals=N (DEFAULT +Inf)  Filter out sites where more than N individuals are heterozygous\n"
"\n"
"       The filtereing parameters that are hard coded are:\n"
"       MIN_OVERALL_VARIANT_PHRED_QUAL=30   The minimum accepted phred scaled [10*log_10*p(variant)] variant quality\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_HELP = 1, OPT_MAX_NUM_HET };

static const char* shortopts = "d:c:s:";

static const struct option longopts[] = {
    { "overall-max-depth",       required_argument, NULL, 'd' },
    { "min-copies",       required_argument, NULL, 'c' },
    { "min-depth-per-sample",       required_argument, NULL, 's' },
    { "max-het-individuals", required_argument, NULL, OPT_MAX_NUM_HET }, 
    { "help",   no_argument, NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int min_copies=1;
    static int max_overall_depth = std::numeric_limits<int>::max();
    static int max_het_indiv = std::numeric_limits<int>::max();
    static int min_depth_in_any_individual = 3;
    static string vcfFile;
}



int testMain(int argc, char** argv) {
    parseTestOptions(argc, argv);
    std::ios_base::openmode mode = std::ios_base::in;
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    std::vector<int> depthsHetFailed;
    std::vector<int> depthsHetPassed;
    std::vector<int> numVariantsPerHetCount;
    std::vector<std::vector<int> > num_indiv_het_vs_depth;
    
    
    std::cerr << "Filtering variants from: " << fileName << std::endl;
    std::cerr << "Maximum read depth (overall) set to: " << opt::max_overall_depth << std::endl;
    std::cerr << "Minimum copies for a variant (e.g. 1 for allowing singletons): " << opt::min_copies << std::endl;
    std::cerr << "Minimum read depth at the variant position: " << opt::min_depth_in_any_individual << std::endl;
    std::cerr << "Filter out sites where more than " << opt::max_het_indiv << " individuals are heterozygous"  << std::endl;
    
    // Open connection to read from the vcf file
    std::ifstream* inFile = new std::ifstream(fileName.c_str(), mode);
    
    
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
                numVariantsPerHetCount.assign(numSamples + 1, 0);
                numChromosomes = (int)numSamples * 2;
                std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
                gotChromosomeNumber = true;
            }
            result.overallQuality = atoi(fields[5].c_str());
            
            
            if (result.overallQuality >= MIN_OVERALL_VARIANT_PHRED_QUAL) {
                result.maxDepthPassed = testOverallReadDepth(opt::max_overall_depth,fields[7]); 
                result.biallelicPassed = testBiallelic(fields[4]);
            }
                
            if (result.biallelicPassed) 
                result.counts = getThisVariantCounts(fields);
            
            
            if (result.biallelicPassed && result.counts.overall <= (numChromosomes - opt::min_copies) && result.counts.overall >= opt::min_copies) {
                if (result.counts.minimumDepthInAnIndividual >= opt::min_depth_in_any_individual) {
                    // filter out sites where more than MAX_NUM_HET individuals are heterozygous
                    bool mnhPassed = testMaxNumHet(result, depthsHetFailed, depthsHetPassed, numVariantsPerHetCount, opt::max_het_indiv,num_indiv_het_vs_depth);
                    if (result.maxDepthPassed && mnhPassed) {
                        std::cout << line << std::endl;
                    }
                }
            }  
        }
    }
    
    std::cerr << "Finished filtering; tabulating depths\n" << std::endl;
    // Does the same as R function table
    std::map<int, int> tableHetFailed = tabulateVector(depthsHetFailed);
    std::map<int, int> tableHetPassed = tabulateVector(depthsHetPassed);
    
    std::ios_base::openmode mode_out = std::ios_base::out;
    string hetFailedFileName = fileRoot + ".het_filter.failed_max" + numToString<int>(opt::max_het_indiv);
    std::ofstream* pHetFailOutFile = new std::ofstream(hetFailedFileName.c_str(), mode_out);
    *pHetFailOutFile << "# Failed HetFilter:" << fileName << std::endl;
    *pHetFailOutFile << "Depth\tNum_variants" << std::endl;
    for (std::map<int, int>::const_iterator iter = tableHetFailed.begin(); iter != tableHetFailed.end(); ++iter) {
        *pHetFailOutFile << iter->first << "\t" << iter->second << std::endl;
    }


    string hetPassedFileName = fileRoot + ".het_filter.passed_max" + numToString<int>(opt::max_het_indiv);
    std::ofstream* pHetPassOutFile = new std::ofstream(hetPassedFileName.c_str(), mode_out);
    *pHetPassOutFile << "# Passed HetFilter:" << fileName << std::endl;
    *pHetPassOutFile << "Depth\tNum_variants" << std::endl;
    for (std::map<int, int>::const_iterator iter = tableHetPassed.begin(); iter != tableHetPassed.end(); ++iter) {
        *pHetPassOutFile << iter->first << "\t" << iter->second << std::endl;
    }
    
    string hetVarCountName = fileRoot + ".het_filter.variants_per_het_count_max" + numToString<int>(opt::max_het_indiv);;
    std::ofstream* pHetVarCountFile = new std::ofstream(hetVarCountName.c_str(), mode_out);
    *pHetVarCountFile << "# Number of variants per each het count: " << fileName << std::endl;
    *pHetVarCountFile << "Num_hets\tNum_variants" << std::endl;
    for (std::vector<int>::size_type i = 0; i != numVariantsPerHetCount.size(); i++) {
        *pHetVarCountFile << i << "\t" << numVariantsPerHetCount[i] << std::endl;
    }
    
    string hetDepthScatterplotName = fileRoot + ".het_filter.depth_scatterplot" + numToString<int>(opt::max_het_indiv);;
    std::ofstream* pHetDepthScatterplotFile = new std::ofstream(hetDepthScatterplotName.c_str(), mode_out);
    *pHetDepthScatterplotFile << "# Number of hets and depth for individual variants for a scatterplot: " << fileName << std::endl;
    *pHetDepthScatterplotFile << "Num_hets\tDepth" << std::endl;
    print_matrix<std::vector<std::vector<int> >&>(num_indiv_het_vs_depth, *pHetDepthScatterplotFile);
    
    return 0;
}




void parseTestOptions(int argc, char** argv) {
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'd': arg >> opt::max_overall_depth; break;
            case 'c': arg >> opt::min_copies; break;
            case 's': arg >> opt::min_depth_in_any_individual; break;
            case '?': die = true; break;
            case OPT_MAX_NUM_HET: arg >> opt::max_het_indiv; break;    
            case OPT_HELP:
                std::cout << TESTING_USAGE_MESSAGE;
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
        std::cout << "\n" << TESTING_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    
}
