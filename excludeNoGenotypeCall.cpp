//
//  excludeNoGenotypeCall.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 28/05/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream> 
#include <map>
#include <getopt.h>
#include "process_vcf_utils.h"

//#define TESTING 1

using std::string;

static const double DIFF_WEIGHT_BOTH_HETS_ME = 0;
static const double DIFF_WEIGHT_BOTH_HETS_RICHARD = (2.0/3.0);
static const double DIFF_WEIGHT_HOM_DIFFERENCE_ME = 1;
static const double DIFF_WEIGHT_HOM_DIFFERENCE_RICHARD = (2.0/3.0);
static const double DIFF_WEIGHT_ONE_HOM_ONE_HET = 0.5;

namespace opt
{
    static int min_copies=1;
    static int max_overall_depth = std::numeric_limits<int>::max();
    static int min_depth_in_any_individual = 3;
    static std::string vcfFile;
    static std::string populationsFile;
    static std::string withBlueIndivFile;
    static bool countHets = false;
    static bool bFilter = true;
    static bool bDiffs = false;
}

enum { OPT_HELP, OPT_DOUBLETON, OPT_HETS, OPT_NOFILTER, OPT_WITH_BLUE, OPT_DIFF_MATRIX };

static const char* shortopts = "d:c:s:";

static const struct option longopts[] = {
    { "overall-max-depth",       required_argument, NULL, 'd' },
    { "min-copies",       required_argument, NULL, 'c' },
    { "min-depth-per-sample",       required_argument, NULL, 's' },
    { "help",   no_argument, NULL, OPT_HELP },
    { "hets-per-individual",   no_argument, NULL, OPT_HETS },
    { "just-stats",   no_argument, NULL, OPT_NOFILTER },
    { "diff-matrix", no_argument,    NULL, OPT_DIFF_MATRIX },
    { "doubleton-distribution", required_argument,    NULL, OPT_DOUBLETON },
    { "count-sites-with-blue", required_argument,    NULL, OPT_WITH_BLUE },
    { NULL, 0, NULL, 0 }
};



// Filtering constants
static const int MIN_OVERALL_VARIANT_PHRED_QUAL=30;   // Corresponds to 0.001 probability of error



static const char *CORRECT_USAGE_MESSAGE =
"Usage: process_vcf [OPTIONS] VCF_FILE\n"
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
"\n"
"       The filtereing parameters that are hard coded are:\n"
"       MIN_OVERALL_VARIANT_PHRED_QUAL=30   The minimum accepted phred scaled [10*log_10*p(variant)] variant quality\n"
"\nReport bugs to mm812@cam.ac.uk\n\n"
"\n"
"       Additional statistics that can be calculated:\n"
"       --doubleton-distribution=POPULATIONS_FILE   Calculates the distribution of doubletons (after all filters) with regards to populations defined in\n"
"                                                   POPULATIONS_FILE; the output matrix is printed to VCF_FILE.doubletons.txt\n"
"       --hets-per-individual                       Outputs the number of heterozygous calls for each individual after passing all\n"
"                                                   filters\n"
"       --just-stats                                Do not do any filtering, just calculate statistics (use in conjunction with the above options)\n"
"                                                   Useful when you want to calculate statistics on an already filtered file\n"
"\nReport bugs to mm812@cam.ac.uk\n\n";



void parseOptions(int argc, char** argv) {
    
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
            case OPT_DIFF_MATRIX: opt::bDiffs = true; break;
            case OPT_DOUBLETON: arg >> opt::populationsFile; break;
            case OPT_HETS: opt::countHets = true; break;
            case OPT_NOFILTER: opt::bFilter = false; break;
            case OPT_WITH_BLUE: arg >> opt::withBlueIndivFile; break;    
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
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
    
    if (!opt::populationsFile.empty() && opt::min_copies > 2) {
        std::cerr << "Can't count doubletons if min-copies > 2\n";
        std::cerr << "Incompatible options\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    
}









int main(int argc, char **argv) {
    typedef std::vector<std::vector<int> >::size_type vec_sz;
    std::ios_base::openmode mode = std::ios_base::in;
    
    if (argc == 1) {
        std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    } else {
        parseOptions(argc, argv);
    }
    
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    std::vector<int> hetCounts;
    bool bDoubletons = !opt::populationsFile.empty();
    bool bCountWithBlue = !opt::withBlueIndivFile.empty();
    std::vector<std::vector<int> > doubletons;
    std::vector<std::vector<double> > diffMatrix;
    std::vector<std::vector<double> > diffMatrixMe;
    std::vector<std::vector<double> > diffMatrixHetsVsHomDiff;
    std::vector<std::string> populations;
    std::vector<std::string> pop_unique;
    std::map<std::string,int> fieldsPopMap;
    MassokoMalawiResult sharingResult;
    
    if (opt::bFilter) {
        std::cerr << "Filtering variants from: " << fileName << std::endl;
        std::cerr << "Maximum read depth (overall) set to: " << opt::max_overall_depth << std::endl;
        std::cerr << "Minimum copies for a variant (e.g. 1 for allowing singletons): " << opt::min_copies << std::endl;
        std::cerr << "Minimum read depth at the variant position: " << opt::min_depth_in_any_individual << std::endl;
    } else {
        std::cerr << "Calculating statistics from: " << fileName << std::endl;
    }    
    if (opt::countHets) {
        
        std::cerr << "Het count per individual are going to be output to: " << fileRoot + ".hets.txt" << std::endl;
    }
    

    std::ifstream* inFile = new std::ifstream(fileName.c_str(), mode);
    
    
    bool gotChromosomeNumber = false;
    int numChromosomes;
    string line;
    int totalVariantNumber = 0;
    while (getline(*inFile, line)) {
        if (line[0] == '#') { if (opt::bFilter) { std::cout << line << std::endl; } } 
        else {
            totalVariantNumber++;
            FilterResult result;
            std::vector<std::string> fields = split(line, '\t');
            if (!gotChromosomeNumber) {
                const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
                numChromosomes = numSamples * 2;
                std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
                gotChromosomeNumber = true;
                hetCounts.assign(numSamples, 0);
            }
            result.overallQuality = atoi(fields[5].c_str());
            
            if (opt::bFilter) {
                if (result.overallQuality >= MIN_OVERALL_VARIANT_PHRED_QUAL) 
                    result.maxDepthPassed = testOverallReadDepth(opt::max_overall_depth,fields[7]); 
                
                if (result.maxDepthPassed)
                    result.biallelicPassed = testBiallelic(fields[4]);
                
                if (result.biallelicPassed) 
                    result.counts = getThisVariantCounts(fields);
                
                if (result.counts.overall <= (numChromosomes - opt::min_copies) && result.counts.overall >= opt::min_copies) {
                    if (result.counts.minimumDepthInAnIndividual >= opt::min_depth_in_any_individual) {
                        std::cout << line << std::endl;
                    }
                }
            }   
        }
    }
    
    
    if (bDoubletons) {
        std::ios_base::openmode mode_out = std::ios_base::out;
        string doubletonFileName = fileRoot + ".doubletons.txt";
        std::ofstream* pDoubletonOutFile = new std::ofstream(doubletonFileName.c_str(), mode_out);
        *pDoubletonOutFile << "# Filters applied:" << std::endl;
        *pDoubletonOutFile << "# --overall-max-depth: " << opt::max_overall_depth << std::endl;
        *pDoubletonOutFile << "# --min-depth-per-sample: " << opt::min_depth_in_any_individual << std::endl;
        *pDoubletonOutFile << "# MIN_OVERALL_VARIANT_PHRED_QUAL=30" << std::endl;
        for (int i = 0; i < pop_unique.size(); i++) {
            if (i == (pop_unique.size()-1))
                *pDoubletonOutFile << pop_unique[i] << std::endl;
            else 
                *pDoubletonOutFile << pop_unique[i] << "\t";
        }
        rearrange_doubletons(doubletons);
        for (int i = 0; i < doubletons.size(); i++) {
            for (int j = 0; j < doubletons.size(); j++) {
                if (j == (doubletons.size()-1))
                    *pDoubletonOutFile << doubletons[i][j] << std::endl;
                else 
                    *pDoubletonOutFile << doubletons[i][j] << "\t";
            }    
        }
    }
    
    if (opt::countHets) {
        std::ios_base::openmode mode_out = std::ios_base::out;
        string hetFileName = fileRoot + ".hets.txt";
        std::ofstream* pHetsOutFile = new std::ofstream(hetFileName.c_str(), mode_out);
        *pHetsOutFile << "# Filters applied:" << std::endl;
        *pHetsOutFile << "# --overall-max-depth: " << opt::max_overall_depth << std::endl;
        *pHetsOutFile << "# --min-copies: " << opt::min_copies << std::endl;
        *pHetsOutFile << "# --min-depth-per-sample: " << opt::min_depth_in_any_individual << std::endl;
        *pHetsOutFile << "# MIN_OVERALL_VARIANT_PHRED_QUAL=30" << std::endl;
        if (bDoubletons) {
            for (int i = 0; i < populations.size(); i++) {
                if (i == (populations.size()-1))
                    *pHetsOutFile << populations[i] << std::endl;
                else 
                    *pHetsOutFile << populations[i] << "\t";
            }
        }
        
        for (int i = 0; i < hetCounts.size(); i++) {
            if (i == (hetCounts.size()-1))
                *pHetsOutFile << hetCounts[i] << std::endl;
            else 
                *pHetsOutFile << hetCounts[i] << "\t";
        }  
    }
    
    if (opt::bDiffs) {
        finalize_diffs_Hets_vs_Homs_proportions(diffMatrixHetsVsHomDiff);
        std::ios_base::openmode mode_out = std::ios_base::out;
        string diffFileName = fileRoot + ".diff_matrix.txt";
        string diffMeFileName = fileRoot + ".diff_me_matrix.txt";
        string hetHomFileName = fileRoot + ".hets_over_homs_matrix.txt";
        std::ofstream* pDiffOutFile = new std::ofstream(diffFileName.c_str(), mode_out);
        std::ofstream* pDiffMeOutFile = new std::ofstream(diffMeFileName.c_str(), mode_out);
        std::ofstream* pHetHomOutFile = new std::ofstream(hetHomFileName.c_str(), mode_out);
        *pDiffOutFile << "# Total number of segragating variant sites in this sample:" << totalVariantNumber << std::endl;
        *pDiffOutFile << "# Richard's scoring scheme" << std::endl;
        *pDiffMeOutFile << "# Total number of segragating variant sites in this sample:" << totalVariantNumber << std::endl;
        *pDiffMeOutFile << "# Homozygous difference = 1, one homozygous, another heterozygous = 0.5:" << totalVariantNumber << std::endl;
        *pHetHomOutFile << "# number of sites both individuals hets/number of sites individuals have a homozygous difference; i.e. num(1/0::1/0)/num(1/1::0/0)" << std::endl;
        *pHetHomOutFile << "# For a free mixing population, we expect this number ~2; for fully separated species ~0" << std::endl;
        if (bDoubletons) {
            for (int i = 0; i < populations.size(); i++) {
                if (i == (populations.size()-1)) {
                    *pDiffOutFile << populations[i] << std::endl;
                    *pDiffMeOutFile << populations[i] << std::endl;
                    *pHetHomOutFile << populations[i] << std::endl;
                } else { 
                    *pDiffOutFile << populations[i] << "\t";
                    *pDiffMeOutFile << populations[i] << "\t";
                    *pHetHomOutFile << populations[i] << "\t";
                }    
            }
        }
        
        for (int i = 0; i < diffMatrix.size(); i++) {
            for (int j = 0; j < diffMatrix.size(); j++) {
                if (j == (diffMatrix.size()-1)) {
                    *pDiffOutFile << diffMatrix[i][j] << std::endl;
                    *pDiffMeOutFile << diffMatrixMe[i][j] << std::endl;
                    *pHetHomOutFile << diffMatrixHetsVsHomDiff[i][j] << std::endl;
                } else {
                    *pDiffOutFile << diffMatrix[i][j] << "\t";
                    *pDiffMeOutFile << diffMatrixMe[i][j] << "\t";
                    *pHetHomOutFile << diffMatrixHetsVsHomDiff[i][j] << "\t";
                }
            }
        }
    }
    
    if (bCountWithBlue) {
        std::cerr << "Getting ready to print" << std::endl;
        std::ios_base::openmode mode_out = std::ios_base::out;
        string shareFileName = fileRoot + ".sharing.txt";
        std::ofstream* pShareOutFile = new std::ofstream(shareFileName.c_str(), mode_out);
        *pShareOutFile << "\t";
        for (int i = 0; i < populations.size(); i++) {
            if (i == (populations.size()-1))
                *pShareOutFile << populations[i] << std::endl;
            else 
                *pShareOutFile << populations[i] << "\t";
        }
        *pShareOutFile << "MassokoFixedBlue(hom/het/absent):\t";
        for (int i = 0; i < populations.size(); i++) {
            if (i == (populations.size()-1))
                *pShareOutFile << sharingResult.homsWithBlue[i] << "/" << sharingResult.hetsWithBlue[i] << "/" << sharingResult.absentWithBlue[i] << std::endl;
            else 
                *pShareOutFile << sharingResult.homsWithBlue[i] << "/" << sharingResult.hetsWithBlue[i] << "/" << sharingResult.absentWithBlue[i] << "\t";
        }
        *pShareOutFile << "MassokoFixedYellow(hom/het/absent):\t";
        for (int i = 0; i < populations.size(); i++) {
            if (i == (populations.size()-1))
                *pShareOutFile << sharingResult.homsWithYellow[i] << "/" << sharingResult.hetsWithYellow[i] << "/" << sharingResult.absentWithYellow[i] << std::endl;
            else 
                *pShareOutFile << sharingResult.homsWithYellow[i] << "/" << sharingResult.hetsWithYellow[i] << "/" << sharingResult.absentWithYellow[i] << "\t";
        }
        
    }
}
