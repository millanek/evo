//
//  evo_getInformativePairs.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 05.04.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#include "evo_getInformativePairs.h"
#include <unordered_set>

#define SUBPROGRAM "InfoPairs"

#define DEBUG 1

static const char *INFOPAIRS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] hapcutBlockFile.txt PAIRTOOLS_FILE.pairs\n"
"Select reads/fragments from the PAIRTOOLS_FILE which could be informative about recombination:\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:";

//enum { OPT_ANNOT, OPT_AF  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string hapcutFile;
    static string pairtoolsFile;
    static string runName = "";
}




int InfoPairsMain(int argc, char** argv) {
    parseInfoPairsOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* hapcutFile = createReader(opt::hapcutFile.c_str());
    std::ifstream* pairtoolsFile = new std::ifstream(opt::pairtoolsFile.c_str());
    
    std::unordered_set<int> phasedHetsPos;
    
    
    // Parse the Hapcut blocks file
    while (getline(*hapcutFile, line)) {
        if (line == "********") continue;
        if (line[0] == 'B' && line[1] == 'L') { // New block - should in the future separate the hets by blocks
            
        } else {
            std::vector<string> phasedSNPdetails = split(line, '\t');
            int snpPos = atoi(phasedSNPdetails[4].c_str());
            phasedHetsPos.insert(snpPos);
        }
    }
    
    // Now parse the pairtools file to find read pairs that can be informative about the phasing and recombination
    while (getline(*pairtoolsFile,line)) {
        // std::cerr << line << std::endl;
        std::vector<string> pairVec = split(line, '\t'); assert(pairVec.size() == 8);
        
        int pair1pos = atoi(pairVec[2].c_str());
        int pair2pos = atoi(pairVec[4].c_str());
        string pair1strand = pairVec[5];
        string pair2strand = pairVec[6];
        
        std::vector<int> posCovered1; posCovered1.resize(150);
        if (pair1strand == "+") std::iota(posCovered1.begin(), posCovered1.end(), pair1pos);
        else std::iota(posCovered1.begin(), posCovered1.end(), pair1pos-150);
        
        std::vector<int> posCovered2; posCovered2.resize(150);
        if (pair1strand == "+") std::iota(posCovered2.begin(), posCovered2.end(), pair2pos);
        else std::iota(posCovered2.begin(), posCovered2.end(), pair2pos-150);
        
        int phasedHetsOnPair = 0;
        
        for (int i = 0; i < posCovered1.size(); i++) {
            if (phasedHetsPos.count(posCovered1[i]) == 1) phasedHetsOnPair++;
        }
        
        for (int i = 0; i < posCovered2.size(); i++) {
            if (phasedHetsPos.count(posCovered2[i]) == 1) phasedHetsOnPair++;
        }
        
        if (phasedHetsOnPair > 1) {
            std::cout << line << std::endl;
        }
        
    }
    
    return 0;
    
}



void parseInfoPairsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'h':
                std::cout << INFOPAIRS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << INFOPAIRS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hapcutFile = argv[optind++];
    opt::pairtoolsFile = argv[optind++];
}
