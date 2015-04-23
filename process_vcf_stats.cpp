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
"       To calculate statistics, it is necessary to supply either:\n"
"       --ind-file=INDIVIDUALS_FILE                 A file listing (one per line) the identities of all the individuals in the vcf file\n"
"       --pop-file=POPULATIONS_FILE                 Defines populations as collections of individuals for statistics that can be calculated per population\n"
"       Then the statistics that can be calculated are:\n"
"       --doubleton-distribution                    Calculates the distribution of doubletons (after all filters) with regards to populations defined in\n"
"                                                   POPULATIONS_FILE; the output matrix is printed to VCF_FILE.doubletons.txt\n"
"       --hets-per-individual                       Outputs the number of heterozygous calls for each individual after passing all\n"
"                                                   filters\n"
"       --diff-matrix                               Generate half-matrices measuring pairwise differences between individuals\n"
"                                                   Useful when you want to calculate statistics on an already filtered file\n"
"\nReport bugs to mm812@cam.ac.uk\n\n";

enum { OPT_INDIV, OPT_POP, OPT_DOUBLETON, OPT_HETS, OPT_DIFF_MATRIX };

static const char* shortopts = "h";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "ind-file", required_argument,    NULL, OPT_INDIV },
    { "pop-file", required_argument,    NULL, OPT_POP },
    { "hets-per-individual",   no_argument, NULL, OPT_HETS },
    { "diff-matrix", no_argument, NULL, OPT_DIFF_MATRIX },
    { "doubleton-distribution", no_argument,    NULL, OPT_DOUBLETON },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    // Individual or population files required for outputtng statistics
    static string individualsFile;
    static string populationsFile;

    // Boolean flags indicating which statistics should be calculated
    static bool countHets = false;
    static bool bDiffs = false;
    static bool bDoubleton = false;
}



int statsMain(int argc, char** argv) {
    parseStatsOptions(argc, argv);
    std::ios_base::openmode mode = std::ios_base::in;
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    
    
    // Data structures to hold individual or population identifiers (headers in the output file(s))
    std::vector<std::string> indPopVector;
    std::vector<std::string> pop_unique;
    std::map<std::string,int> fieldsPopMap;
    // Read in INDIVIDUALS_FILE or POPULATIONS_FILE
    std::ifstream* indOrPopFile = !opt::individualsFile.empty() ? new std::ifstream(opt::individualsFile.c_str(), mode) : new std::ifstream(opt::populationsFile.c_str(), mode);
    string indPopLine; 
    while (getline(*indOrPopFile, indPopLine)) {
        indPopVector.push_back(indPopLine);
    }
    
    // Data structures to hold the results
    std::vector<int> hetCounts;
    std::vector<std::vector<int> > doubletons;
    std::vector<std::vector<double> > diffMatrix;
    std::vector<std::vector<double> > diffMatrixMe;
    std::vector<std::vector<double> > diffMatrixHetsVsHomDiff;
    initialize_matrix_double(diffMatrix, (int)indPopVector.size());
    initialize_matrix_double(diffMatrixMe, (int)indPopVector.size());
    initialize_matrix_double(diffMatrixHetsVsHomDiff, (int)indPopVector.size());
    
    
    std::cerr << "Calculating statistics from: " << fileName << std::endl;
    if (opt::countHets) 
        std::cerr << "Het count per individual are going to be output to: " << fileRoot + ".hets.txt" << std::endl;
    if (opt::bDoubleton) { 
        std::cerr << "The distribution of doubletons is going to be output to: " <<  fileRoot + ".doubletons.txt" << std::endl;
        pop_unique = initializeDoubletons(doubletons, indPopVector,fieldsPopMap);
    }
    
    // Start reading from the vcf file
    std::istream* inFile = createReader(fileName.c_str());
    bool gotChromosomeNumber = false;
    int numChromosomes;
    string line;
    int totalVariantNumber = 0;
    while (getline(*inFile, line)) {
        if (line[0] != '#') {  
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
            
            result.counts = getThisVariantCounts(fields);
            if (opt::countHets) {
                het_analysis(hetCounts,result);
            }
            if (opt::bDoubleton) {
                doubleton_analysis(doubletons,result,numChromosomes,indPopVector, fieldsPopMap);  
            }
            if (opt::bDiffs) {
                diffs_between_individuals(diffMatrix,diffMatrixMe,diffMatrixHetsVsHomDiff,result);
            }
        }
    }
    
    // Printing doubletons
    if (opt::bDoubleton) 
        print_doubleton_distribution(fileRoot, pop_unique, doubletons);
    
    // Printing het counts
    if (opt::countHets) 
        print_het_counts(fileRoot, indPopVector, hetCounts);  
    
    // Printing pairwise difference statistics
    if (opt::bDiffs) {
        finalize_diffs_Hets_vs_Homs_proportions(diffMatrixHetsVsHomDiff);
        print_pairwise_diff_stats(fileRoot, indPopVector, totalVariantNumber, diffMatrix, diffMatrixMe, diffMatrixHetsVsHomDiff);
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
            case OPT_INDIV: arg >> opt::individualsFile; break;    
            case OPT_DIFF_MATRIX: opt::bDiffs = true; break;    
            case OPT_DOUBLETON: opt::bDoubleton = true; break;
            case OPT_HETS: opt::countHets = true; break; 
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
    
    // Make sure we have INDIVIDUALS_FILE or POPULATIONS_FILE if calculating statistics
    if (opt::populationsFile.empty() && opt::individualsFile.empty()) {
        std::cerr << "Cannot calculate statistics without INDIVIDUALS_FILE or POPULATIONS_FILE\n";
        die = true;
    }
    
    if (!opt::bDiffs && !opt::bDoubleton && !opt::countHets) {
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


