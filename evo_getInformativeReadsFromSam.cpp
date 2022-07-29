//
//  evo_getInformativeReadsFromSam.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 05.04.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//


// Sam flags:
// 65 -     Paired, first in pair, + strand
// 81 -     Paired, first in pair, - strand
// 97 -     Paired, first in pair, + strand
// 113 -    Paired, first in pair, - strand
// 129 - Paired, second in pair, + strand
// 145 - Paired, second in pair, - strand
// 161 - Paired, second in pair, + strand
// 177 - Paired, second in pair, - strand
// >2000 - secondary alignment


#include "evo_getInformativeReadsFromSam.h"
#include <unordered_set>

#define SUBPROGRAM "InfoReadsSam"

#define DEBUG 1

static const char *INFOREADS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] hetPosFile.txt\n"
"Select reads from SAMTOOLS_FILE which could be informative about recombination\n"
"Expects sam input on STDIN\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -m, --min-MQ                            (default: 20) the minimum mapping quality for a read to be considered\n"
"       --hapCut                                the het positions come from HapCut output"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:m:";

enum { OPT_HAPCUT  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "min-MQ",   required_argument, NULL, 'm' },
    { "run-name",   required_argument, NULL, 'n' },
    { "hapCut",   no_argument, NULL, OPT_HAPCUT },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string hetFile;
    static bool hapcutFormat = false;
    static string runName = "";
    static int minMQ = 20;
}




int InfoReadsMain(int argc, char** argv) {
    parseInfoReadsOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* hetFile = createReader(opt::hetFile.c_str());
    std::istream* samtoolsFile = &std::cin;
    
    std::unordered_set<int> phasedHetsPos;
    std::map<int,PhaseInfo*> positionToPhase;
    
    
    if (opt::hapcutFormat) {
        // Parse the Hapcut blocks file
        while (getline(*hetFile, line)) {
            if (line[0] == '*') {
            
            } else if (line[0] == 'B' && line[1] == 'L') { // New block - should in the future separate the hets by blocks
                
            } else {
                std::vector<string> phasedSNPdetails = split(line, '\t');
                int snpPos = atoi(phasedSNPdetails[4].c_str());
                int H1phase = atoi(phasedSNPdetails[1].c_str());
                int H2phase = atoi(phasedSNPdetails[2].c_str());
                //std::cout << "line: " << line << std::endl;
                // std::cout << "phasedSNPdetails[5]: " << phasedSNPdetails[5] << std::endl;
                assert(phasedSNPdetails[5].length() == 1); assert(phasedSNPdetails[6].length() == 1);
                char refBase = phasedSNPdetails[5][0];
                char altBase = phasedSNPdetails[6][0];
                double phaseQual = stringToDouble(phasedSNPdetails[10].c_str());
                int snpCoverage = atoi(phasedSNPdetails[11].c_str());
                std::vector<char> phasedVars;
                if (H1phase == 0 && H2phase == 1) {
                    phasedVars.push_back(refBase); phasedVars.push_back(altBase);
                } else if (H1phase == 1 && H2phase == 0) {
                    phasedVars.push_back(altBase); phasedVars.push_back(refBase);
                } else{
                    continue;
                }
                PhaseInfo* thisPhase = new PhaseInfo(snpPos,phaseQual,snpCoverage, phasedVars);
                positionToPhase[snpPos] = thisPhase;
            }
        }
    } else {
        while (getline(*hetFile, line)) {
            std::vector<string> phasedSNPdetails = split(line, '\t');
            int snpPos = atoi(phasedSNPdetails[1].c_str());
            phasedHetsPos.insert(snpPos);
        }
    }
    
    // Now parse the pairtools file to find read pairs that can be informative about the phasing and recombination
    while (getline(*samtoolsFile,line)) {
        // std::cerr << line << std::endl;
        std::vector<string> samRecVec = split(line, '\t'); //assert(pairVec.size() == 8);
        int flag = atoi(samRecVec[1].c_str());
        if (flag > 2000) continue;
        
        int MQ = atoi(samRecVec[4].c_str());
        if (MQ < opt::minMQ) continue;
        
        
        RecombRead* thisRead = new RecombRead(samRecVec);
        
        thisRead->findHetsInRead(positionToPhase);
        
        if (thisRead->hetSites.size() > 0) {
            std::cout << line << std::endl;
        }
    }
    
    return 0;
    
}



void parseInfoReadsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case OPT_HAPCUT: opt::hapcutFormat = true; break;
            case 'm': arg >> opt::minMQ; break;
            case 'h':
                std::cout << INFOREADS_USAGE_MESSAGE;
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
        std::cout << "\n" << INFOREADS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetFile = argv[optind++];
}
