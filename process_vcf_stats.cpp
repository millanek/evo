//
//  process_vcf_stats.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_stats.h"
#include "process_vcf_stats_functions.h"
#include "process_vcf_print_routines.h"

#define SUBPROGRAM "stats"

static const char *STATS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE\n"
"Calculating various statistics from a VCF file, output to STD_OUT\n"
"\n"
"       --help                           display this help and exit\n"
"\n"
"       --ind-file=INDIVIDUALS_FILE                 A file listing (one per line) the identities of all the individuals in the vcf file\n"
"                                                   (if not supplied, individual IDs are read from the VCF header)\n"
"       To calculate statistics for populations, it is necessary to supply:\n"
"       --pop-file=POPULATIONS_FILE                 Defines populations as collections of individuals for statistics that can be calculated per population\n"
"       Then the statistics that can be calculated are:\n"
"       --doubleton-distribution                    Calculates the distribution of doubletons (after all filters) with regards to populations defined in\n"
"       --private-variants                          Calculates the numbers of private variants in each of the populations defined in\n"
"                                                   POPULATIONS_FILE; the output matrix is printed to VCF_FILE.privateVars.txt\n"
"       --hets-per-individual                       Outputs the number of heterozygous calls for each individual after passing all\n"
"                                                   filters\n"
"       --diff-matrix                               Generate half-matrices measuring pairwise differences between individuals\n"
"                                                   Useful when you want to calculate statistics on an already filtered file\n"
"       --diff-matrix-h1                            Generate a half-matrix measuring pairwise differences between haplotypes1\n"
"       --diff-matrix-allH                          Generate a half-matrix measuring pairwise differences between all haplotypes\n"
"       The program can also output 100 bootstrap replicates of the distance matrices unsing a block 'case resampling' scheme:\n"
"       --block-bootstrap=BLOCKSIZE (default 100)  Generate 100 distance matrices resampling with replacement in blocks of BLOCKSIZE variants\n"
"\nReport bugs to" PACKAGE_BUGREPORT "\n\n";

enum { OPT_INDIV, OPT_POP, OPT_DOUBLETON, OPT_HETS, OPT_DIFF_MATRIX, OPT_DIFF_MATRIX_H1, OPT_DIFF_MATRIX_ALLH, OPT_BLOCK_BOOTSTRAP, OPT_PRIVATE_VARS };

static const char* shortopts = "h";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "ind-file", required_argument,    NULL, OPT_INDIV },
    { "pop-file", required_argument,    NULL, OPT_POP },
    { "hets-per-individual",   no_argument, NULL, OPT_HETS },
    { "diff-matrix", no_argument, NULL, OPT_DIFF_MATRIX },
    { "doubleton-distribution", no_argument,    NULL, OPT_DOUBLETON },
    { "diff-matrix-h1", no_argument,    NULL, OPT_DIFF_MATRIX_H1 },
    { "diff-matrix-allH", no_argument,    NULL, OPT_DIFF_MATRIX_ALLH },
    { "block-bootstrap", required_argument,    NULL, OPT_BLOCK_BOOTSTRAP },
    { "private-variants", required_argument,    NULL, OPT_PRIVATE_VARS },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    // Individual or population files required for outputtng statistics
    static string sampleNameFile;
    static string populationsFile = "";

    // Boolean flags indicating which statistics should be calculated
    static bool countHets = false;
    static bool countPrivateVars = false;
    static bool bDiffs = false;
    static bool bDoubleton = false;
    static bool bDiffH1 = false;
    static bool bDiffAllH = false;
    
    static int bootstrapBlockSize = 0;
}

int statsMain(int argc, char** argv) {
    parseStatsOptions(argc, argv);
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    std::ifstream* popFile;
    
    // Data structures to hold individual and population identifiers
    std::vector<std::string> sampleNames; std::vector<string> populationLabels;
    std::vector<std::string> populationsStrings; std::vector<std::vector<size_t> > populationsIndices;
    // Read in the POPULATIONS_FILE if supplied
    if (!opt::populationsFile.empty()) {
        popFile = new std::ifstream(opt::populationsFile.c_str());
        string line;
        while (getline(*popFile, line)) {
            std::vector<std::string> popnamePopstring = split(line, '\t');
            populationLabels.push_back(popnamePopstring[0]);
            populationsStrings.push_back(popnamePopstring[1]);
        }
    }
    
    // Data structures to hold the results
    std::vector<int> hetCounts; std::vector<int> hetsSharedWithOthers; std::vector<int> privateVarCounts;
    std::vector<std::vector<int> > doubletons; std::vector<std::vector<double> > diffMatrixHetsVsHomDiff;
    std::vector<std::vector<double> > diffMatrix; std::vector<std::vector<double> > diffMatrixMe;
    std::vector<std::vector<double> > diffMatrixH1; std::vector<std::vector<double> > diffMatrixAllH;
    
    std::cerr << "Calculating statistics from: " << fileName << std::endl;
    if (opt::countHets) {
        std::cerr << "Het count per individual are going to be output to: " << fileRoot + ".hets.txt" << std::endl;
        std::cerr << "Counts of hets (per individual) shared with others are going to be output to: " << fileRoot + ".sharedHets.txt" << std::endl;
    }
    if (opt::bDoubleton) {
        std::cerr << "The distribution of doubletons is going to be output to: " <<  fileRoot + ".doubletons.txt" << std::endl;
        //pop_unique = initializeDoubletons(doubletons, indPopVector,fieldsPopMap);
    }
    if (opt::countPrivateVars) {
        std::cerr << "The numbers of private variants fixed in groups defined in " << opt::populationsFile << " are going to be output to: " <<  fileRoot + ".privateFixedVars.txt" << std::endl;
        //pop_unique = initializeDoubletons(doubletons, indPopVector,fieldsPopMap);
    }
    
    // Start reading from the vcf file
    std::istream* inFile = createReader(fileName.c_str());
    int numSamples; int numChromosomes; int totalVariantNumber = 0;
    string line;
    while (getline(*inFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            // Read the sample names, initialise output variables, and (if supplied) locate populations
            std::vector<std::string> fields = split(line, '\t');
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            } numSamples = (int)sampleNames.size(); numChromosomes = (int)numSamples * 2;
            std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
            initialize_matrix_double(diffMatrix, numSamples); initialize_matrix_double(diffMatrixMe, numSamples);
            initialize_matrix_double(diffMatrixHetsVsHomDiff, numSamples);
            initialize_matrix_double(diffMatrixH1, numSamples); initialize_matrix_double(diffMatrixAllH, numSamples);
            privateVarCounts.assign(numSamples, 0); hetCounts.assign(numSamples, 0); hetsSharedWithOthers.assign(numSamples, 0);
            if (!populationsStrings.empty()) {
                for (int i = 0; i != (int)populationsStrings.size(); i++) {
                    std::vector<size_t> thisIndices = locateSet(sampleNames, split(populationsStrings[i],','));
                    populationsIndices.push_back(thisIndices);
                }
            }
        } else {
            totalVariantNumber++;
            std::vector<std::string> fields = split(line, '\t');
            
            FilterResult result;
            result.counts = getThisVariantCounts(fields);
            
            // Only do these calculations if none of the genotypes are missing:
            if (result.counts.bAnyMissingGenotypes == false) {
                if (opt::countHets) {
                    het_analysis(hetCounts, hetsSharedWithOthers, result);
                }
                if (opt::countPrivateVars) {
                    privateVars_analysis(privateVarCounts,result,populationsIndices);
                }
                if (opt::bDoubleton) {
                //    doubleton_analysis(doubletons,result,numChromosomes,indPopVector, fieldsPopMap);
                }
                if (opt::bDiffs) {
                    diffs_between_individuals(diffMatrix,diffMatrixMe,diffMatrixHetsVsHomDiff,result);
                }
                if (opt::bDiffH1) {
                    diffs_between_H1(diffMatrixH1, result);
                }
                if (opt::bDiffAllH) {
                    diffs_between_AllH(diffMatrixAllH, result);
                }
            }
            if (totalVariantNumber % 100000 == 0)
                std::cerr << "Processed " << totalVariantNumber << " variants" << std::endl;
            
        }
    }
    
    // Printing doubletons
//    if (opt::bDoubleton)
//        print_doubleton_distribution(fileRoot, pop_unique, doubletons);
    
    // Printing het counts
    if (opt::countHets) 
        print_het_counts(fileRoot, sampleNames, hetCounts, hetsSharedWithOthers);
    
    // Printing statistics for private variants
    if (opt::countPrivateVars)
    print_privateFixedVarsSummary(fileRoot, populationLabels, opt::populationsFile, privateVarCounts);
    
    // Printing pairwise difference statistics
    if (opt::bDiffs) {
        finalize_diffs_Hets_vs_Homs_proportions(diffMatrixHetsVsHomDiff);
        print_pairwise_diff_stats(fileRoot, sampleNames, totalVariantNumber, diffMatrix, diffMatrixMe, diffMatrixHetsVsHomDiff);
    }
    if (opt::bDiffH1) {
        print_H1_pairwise_diff_stats(fileRoot, sampleNames, totalVariantNumber, diffMatrixH1);
    }
    if (opt::bDiffAllH) {
        print_AllH_pairwise_diff_stats(fileRoot, sampleNames, totalVariantNumber, diffMatrixAllH);
    }
    
    return 0;
}

void parseStatsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case OPT_POP: arg >> opt::populationsFile; break;
            case OPT_INDIV: arg >> opt::sampleNameFile; break;
            case OPT_DIFF_MATRIX: opt::bDiffs = true; break;    
            case OPT_DOUBLETON: opt::bDoubleton = true; break;
            case OPT_HETS: opt::countHets = true; break;
            case OPT_DIFF_MATRIX_H1: opt::bDiffH1 = true; break;
            case OPT_DIFF_MATRIX_ALLH: opt::bDiffAllH = true; break;
            case OPT_BLOCK_BOOTSTRAP: arg >> opt::bootstrapBlockSize; break;
            case OPT_PRIVATE_VARS: opt::countPrivateVars = true; break;
            case 'h':
                std::cout << STATS_USAGE_MESSAGE;
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
    
    if (opt::countPrivateVars && opt::populationsFile == "") {
        std::cerr << "You need to supply a file defining populations to count private fixed variants\n";
        die = true;
    }
    
    if (!opt::bDiffs && !opt::bDoubleton && !opt::countHets && !opt::countPrivateVars) {
        std::cerr << "Which statistics should the program calculate?\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << STATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
}


