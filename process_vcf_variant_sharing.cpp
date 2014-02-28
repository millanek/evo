//
//  process_vcf_variant_sharing.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_variant_sharing.h"


#define SUBPROGRAM "variant-sharing"

static const char *SHARING_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE\n"
"Comparing variant segragation between Massoko and Malawi, output to STD_OUT\n"
"\n"
"       --help                           display this help and exit\n"
"\n"
"       To calculate statistics, it is necessary to supply either:\n"
"       --ind-file=INDIVIDUALS_FILE                 A file listing (one per line) the identities of all the individuals in the vcf file\n"
"       --pop-file=POPULATIONS_FILE                 Defines populations as collections of individuals for statistics that can be calculated per population\n"
"\nReport bugs to mm812@cam.ac.uk\n\n";

namespace opt
{
    static string vcfFile;
    
    // Individual or population files required for outputtng statistics
    static string individualsFile;
    static string populationsFile;
    static string withBlueIndivFile;
}

enum { OPT_HELP, OPT_WITH_BLUE, OPT_INDIV, OPT_POP };

static const char* shortopts = "h";

static const struct option longopts[] = {
    { "overall-max-depth",       required_argument, NULL, 'd' },
    { "min-copies",       required_argument, NULL, 'c' },
    { "min-depth-per-sample",       required_argument, NULL, 's' },
    { "help",   no_argument, NULL, 'h' },
    { "ind-file", required_argument,    NULL, OPT_INDIV },
    { "pop-file", required_argument,    NULL, OPT_POP },
    { "count-sites-with-blue", required_argument,    NULL, OPT_WITH_BLUE },
    { NULL, 0, NULL, 0 }
};


int variantSharingMain(int argc, char** argv) {
    parseVariantSharingOptions(argc, argv);
    std::ios_base::openmode mode = std::ios_base::in;
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    
    // For variant sharing between Massoko and Malawi
    bool bCountWithBlue = !opt::withBlueIndivFile.empty();
    MassokoMalawiResult sharingResult;
    
    // Data structures to hold individual or population identifiers (headers in the output file(s))
    std::vector<std::string> populations;
    std::vector<std::string> pop_unique;
    std::map<std::string,int> fieldsPopMap;
    
    std::cerr << "Analysing sharing of fixed variants between Massoko and Lake Malawi" << std::endl;
    std::cerr << "The ind file should be: sample_individuals_All.txt" << std::endl;
    std::ifstream* indFile = new std::ifstream(opt::withBlueIndivFile.c_str(), mode);
    string ind; 
    while (getline(*indFile, ind)) {
        populations.push_back(ind);
    }
    if (populations.size() != 46) {
        std::cerr << "It seems you supplied wrong file; the correct file should have 46 individuals" << std::endl;
        exit(EXIT_FAILURE);
    }
    populations.erase(populations.begin(), populations.begin()+12); // Remove the Massoko fish from this list of individuals
    sharingResult.init(populations.size());
    
    // Start reading from the vcf file
    std::ifstream* inFile = new std::ifstream(fileName.c_str(), mode);
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
            }
            result.overallQuality = atoi(fields[5].c_str());
            result.counts = getThisVariantCounts(fields);
            if (bCountWithBlue) {
                massokoMalawiSharing(result, sharingResult);
            }   
        }
    }
    
    // Print out the results
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
        
    return 0;
}


void parseVariantSharingOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case OPT_POP: arg >> opt::populationsFile; break;
            case OPT_INDIV: arg >> opt::individualsFile; break;    
            case 'h':
                std::cout << SHARING_USAGE_MESSAGE;
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
        std::cerr << "You should provide an INDIVIDUALS_FILE or POPULATIONS_FILE (most likely sample_individuals_All.txt)\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << SHARING_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
}