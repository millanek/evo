//
//  process_vcf_massoko.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 11/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_massoko.h"

// Filtering constants
static const int MIN_OVERALL_VARIANT_PHRED_QUAL=30;   // Corresponds to 0.001 probability of error

#define SUBPROGRAM "massoko"

static const char *FILTER_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE\n"
"Custom file for filtering and getting fixed (and nearly fixed) variants from the Lake Massoko VCF file, output to STD_OUT\n"
"\n"
"       --help                           display this help and exit\n"
"\n"
"       The filtereing parameters that can be changed are:\n"
"       -d, --overall-max-depth (DEFAULT +Inf)  Maximum read depth allowed at the putative variant site - filters out variants due to\n"
"                                               collapsed repeats in the reference (reads from multiple sites are all mapped to one)\n"
"       -c, --min-copies=MIN (DEFAULT 1)        The variant needs to be present in at least MIN copies\n"
"                                               (i.e. setting --min_copies==1 gives all segregating sites including singletons)\n"
"       -s, --min-depth-per-sample=MIN (DEFAULT 3)        Minimum read depth at the variant position for any sample\n"
"                                               (i.e. all samples neet to have at least MIN reads at the variant position)\n"
"\n"
"       The filtereing parameters that are hard coded are:\n"
"       MIN_OVERALL_VARIANT_PHRED_QUAL=30   The minimum accepted phred scaled [10*log_10*p(variant)] variant quality\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "d:c:s:";

static const struct option longopts[] = {
    { "overall-max-depth",       required_argument, NULL, 'd' },
    { "min-copies",       required_argument, NULL, 'c' },
    { "min-depth-per-sample",       required_argument, NULL, 's' }, 
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

class MassokoCounts : public Counts {
public:
    MassokoCounts() : Counts(), blue(0) {};
    
    int blue;
};

MassokoCounts getMassokoVariantCounts(const std::vector<std::string>& fields) {
    MassokoCounts thisVariantCounts;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') { 
            thisVariantCounts.overall++;
            if (i <= 14) { thisVariantCounts.blue++; }
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') { 
            thisVariantCounts.overall++; 
            if (i <= 14) { thisVariantCounts.blue++; }
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
        std::vector<std::string> genotypeData = split(fields[i], ':');
        
        if (atoi(genotypeData[2].c_str()) < thisVariantCounts.minimumDepthInAnIndividual) {
            thisVariantCounts.minimumDepthInAnIndividual = atoi(genotypeData[2].c_str());
        }
        
    }
    // Also get overall depth for this variant
    std::vector<std::string> info = split(fields[7], ';');
    if (info[0] == "INDEL") {
        info = split(info[1], '=');
    } else {
        info = split(info[0], '=');
    }
    thisVariantCounts.overallDepth = atoi((info.back()).c_str());
    return thisVariantCounts;
}

// Results object
class MassokoFilterResult
{
public:
    MassokoFilterResult() : overallQuality(0), maxDepthPassed(false), counts(), biallelicPassed(false), degenPassed(false) {}
    
    int overallQuality;  // -10log_10 p(no variant)
    bool maxDepthPassed;
    MassokoCounts counts;
    bool biallelicPassed;
    bool degenPassed;
};

void initializeCountsBlue(std::map<int, std::vector<int> >& countsBlue, int numChr) {
    for (int i = opt::min_copies; i <= numChr - opt::min_copies; i++) {
        std::vector<int> v(i+1);
        for (int j = 0; j <= i; j++) { v[j] = 0; }
        countsBlue[i] = v;
    }
}

int getNumHets(MassokoFilterResult& result) {
    int num_hets = 0;
    for (std::vector<std::vector<int> >::size_type i = 0; i < result.counts.individualsWithVariant.size(); i++) {
        if (result.counts.individualsWithVariant[i] == 1)
            num_hets++;
    }
    return num_hets;
}

int massokoMain(int argc, char** argv) {
    parseMassokoOptions(argc, argv);
    std::ios_base::openmode mode = std::ios_base::in;
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    std::map<int, std::vector<int> > counts_blue;
    
    std::cerr << "Filtering variants from: " << fileName << std::endl;
    std::cerr << "Maximum read depth (overall) set to: " << opt::max_overall_depth << std::endl;
    std::cerr << "Minimum copies for a variant (e.g. 1 for allowing singletons): " << opt::min_copies << std::endl;
    std::cerr << "Minimum read depth at the variant position: " << opt::min_depth_in_any_individual << std::endl;
    
    // Open connection to read from the vcf file
    std::ifstream* inFile = new std::ifstream(fileName.c_str(), mode);
    std::ios_base::openmode mode_out = std::ios_base::out;
    string VarFile12FileName = fileRoot + ".fixed_variants12.txt";
    string VarFile11FileName = fileRoot + ".fixed_variants11.txt";
    string VarFile10FileName = fileRoot + ".fixed_variants10_two_hets.txt";
    string VarFile10HomFileName = fileRoot + ".fixed_variants10_hom.txt";
    string VarFile9FileName = fileRoot + ".fixed_variants9_three_hets.txt";
    string VarFile9HomFileName = fileRoot + ".fixed_variants9_hom.txt";
    std::ofstream* pFixedVarFile12 = new std::ofstream(VarFile12FileName.c_str(), mode_out);
    std::ofstream* pFixedVarFile11 = new std::ofstream(VarFile11FileName.c_str(), mode_out);
    std::ofstream* pFixedVarFile10 = new std::ofstream(VarFile10FileName.c_str(), mode_out);
    std::ofstream* pFixedVarFile10hom = new std::ofstream(VarFile10HomFileName.c_str(), mode_out);
    std::ofstream* pFixedVarFile9 = new std::ofstream(VarFile9FileName.c_str(), mode_out);
    std::ofstream* pFixedVarFile9hom = new std::ofstream(VarFile9HomFileName.c_str(), mode_out);
    bool gotChromosomeNumber = false;
    int numChromosomes;
    string line;
    int totalVariantNumber = 0;
    while (getline(*inFile, line)) {
        if (line[0] == '#') {
            std::cout << line << std::endl;
            *pFixedVarFile12 << line << std::endl; 
            *pFixedVarFile11 << line << std::endl;
            *pFixedVarFile10 << line << std::endl;
        } else {
            totalVariantNumber++;
            MassokoFilterResult result;
            std::vector<std::string> fields = split(line, '\t');
            if (!gotChromosomeNumber) {
                const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
                numChromosomes = numSamples * 2;
                std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
                gotChromosomeNumber = true;
                initializeCountsBlue(counts_blue, numChromosomes);
            }
            result.overallQuality = atoi(fields[5].c_str());
            
            
            if (result.overallQuality >= MIN_OVERALL_VARIANT_PHRED_QUAL) 
                result.maxDepthPassed = testOverallReadDepth(opt::max_overall_depth,fields[7]); 
            
            if (result.maxDepthPassed)
                result.biallelicPassed = testBiallelic(fields[4]);
            
            if (result.biallelicPassed) 
                result.counts = getMassokoVariantCounts(fields);
            
            // filter out sites where more than MAX_NUM_HET individuals are heterozygous
            int num_hets = getNumHets(result);
            
            if (result.biallelicPassed && result.counts.overall <= (numChromosomes - opt::min_copies) && result.counts.overall >= opt::min_copies) {
                if (result.counts.minimumDepthInAnIndividual >= opt::min_depth_in_any_individual) {
                    std::cout << line << std::endl;
                    counts_blue[result.counts.overall][result.counts.blue]++;
                    if (result.counts.overall == 12 && (result.counts.blue == 12 || result.counts.blue == 0)) {
                        *pFixedVarFile12 << line << std::endl;
                    }
                    if (result.counts.overall == 11 && (result.counts.blue == 11 || result.counts.blue == 0)) {
                        *pFixedVarFile11 << line << std::endl;
                    }
                    if (result.counts.overall == 13 && (result.counts.blue == 12 || result.counts.blue == 1)) {
                        *pFixedVarFile11 << line << std::endl;
                    }
                    if (result.counts.overall == 10 && (result.counts.blue == 10 || result.counts.blue == 0) && (num_hets == 2)) {
                        *pFixedVarFile10 << line << std::endl;
                    } else if (result.counts.overall == 10 && (result.counts.blue == 10 || result.counts.blue == 0)) {
                        *pFixedVarFile10hom << line << std::endl;
                    }
                    if (result.counts.overall == 14 && (result.counts.blue == 12 || result.counts.blue == 2) && (num_hets == 2)) {
                        *pFixedVarFile10 << line << std::endl;
                    } else if (result.counts.overall == 14 && (result.counts.blue == 14 || result.counts.blue == 2)) {
                        *pFixedVarFile10hom << line << std::endl;
                    }
                    if (result.counts.overall == 9 && (result.counts.blue == 9 || result.counts.blue == 0) && (num_hets == 3)) {
                        *pFixedVarFile9 << line << std::endl;
                    } else if (result.counts.overall == 9 && (result.counts.blue == 9 || result.counts.blue == 0)) {
                        *pFixedVarFile9hom << line << std::endl;
                    }
                    if (result.counts.overall == 15 && (result.counts.blue == 12 || result.counts.blue == 3) && (num_hets == 3)) {
                        *pFixedVarFile9 << line << std::endl;
                    } else if (result.counts.overall == 15 && (result.counts.blue == 12 || result.counts.blue == 3)) {
                        *pFixedVarFile9hom << line << std::endl;
                    }
                }
            }
        }
    }
    for (int i = numChromosomes - opt::min_copies; i >= opt::min_copies; i--) {
        for (std::vector<std::string>::size_type j = 0; j != counts_blue[i].size(); j++) {
            if (j == (counts_blue[i].size()-1)) {
                std::cerr << counts_blue[i][j];
            } else {
                std::cerr << counts_blue[i][j] << "\t";
            }
        } std::cerr << std::endl;
    } 
    return 0;
}

void parseMassokoOptions(int argc, char** argv) {
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

