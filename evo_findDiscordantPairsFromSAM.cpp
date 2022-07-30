//
//  evo_findDiscordantPairs.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 08.04.22.
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


#include "evo_findDiscordantPairsFromSAM.h"

#define SUBPROGRAM "DiscordantPairsFromSAM"

#define DEBUG 1

static const char *DISCORDPAIRS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] hapcutBlockFile.txt INFORMATVE_PAIRS.sam\n"
"Select reads/fragments from the PAIRTOOLS_FILE which could be informative about recombination:\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -m, --min-MQ                            (default: 20) the minimum mapping quality for a read to be considered\n"
"       -b, --min-BQ                            (default: 30) the minimum base quality for assesssing discordant phase\n"
"       -p, --min-PQ                            (default: 30) the minimum phase quality for assesssing discordant phase (relevant with the --hapCut option)\n"
"       --hapCut                                the het positions come from HapCut output\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:b:m:p:";

//enum { OPT_ANNOT, OPT_AF  };
enum { OPT_HAPCUT  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "min-MQ",   required_argument, NULL, 'm' },
    { "min-BQ",   required_argument, NULL, 'b' },
    { "min-PQ",   required_argument, NULL, 'p' },
    { "hapCut",   no_argument, NULL, OPT_HAPCUT },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static bool hapcutFormat = false;
    static string hetsFile;
    static string samFile;
    static string runName = "";
    static int minMQ = 20;
    static int minBQ = 30;
    static int minPQ = 30;
}




int DiscordPairsFromSAMMain(int argc, char** argv) {
    parseDiscordPairsFromSAMOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* hetsFile = createReader(opt::hetsFile.c_str());
    std::ifstream* samFile = new std::ifstream(opt::samFile.c_str());
    
    std::ofstream* phaseSwitchFile = new std::ofstream("switches" + opt::runName + ".txt");
    std::ofstream* goodReadPairsFile = new std::ofstream("goodReadPairs" + opt::runName + ".txt");
    
    std::map<int,PhaseInfo*> positionToPhase;
    std::map<string,std::vector<int>> infoPairNameToPos;
    std::map<string,std::vector<string>> infoPairNameToStrands;
    
    std::map<string,std::vector<RecombRead*>> samNameToReads;
    int numHetPairs = 0;
    std::unordered_map<string, ReadLinkSNPpair*> SNPpairs;
    
    if (opt::hapcutFormat) {
        // Parse the Hapcut blocks file
        while (getline(*hetsFile, line)) {
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
        while (getline(*hetsFile, line)) {
            std::vector<string> phasedSNPdetails = split(line, '\t');
            int snpPos = atoi(phasedSNPdetails[1].c_str());
            assert(phasedSNPdetails[2].length() == 1); assert(phasedSNPdetails[3].length() == 1);
            char refBase = phasedSNPdetails[2][0];
            char altBase = phasedSNPdetails[3][0];
            
            std::vector<char> phasedVars;
            phasedVars.push_back(refBase); phasedVars.push_back(altBase);
            double phaseQual = 0; int snpCoverage = 0;
            PhaseInfo* thisPhase = new PhaseInfo(snpPos,phaseQual,snpCoverage, phasedVars);
            positionToPhase[snpPos] = thisPhase;
        }
        numHetPairs = nChoosek((int)positionToPhase.size(),2);
    }
    
    
    
    
    
    // Now parse the samtools file to find reads that match records from the pairstools file
    // and therefore  can be informative about the phasing and recombination
    int readN = 0;
    std::vector<RecombRead*> informativeReads;
    while (getline(*samFile,line)) {
        // std::cerr << line << std::endl;
        readN++;
        std::vector<string> samRecVec = split(line, '\t'); //assert(pairVec.size() == 8);
        RecombRead* thisRead = new RecombRead(samRecVec);
        informativeReads.push_back(thisRead);
    }
    
    std::vector<RecombReadPair*> informativeReadPairs;
    for (int r = 0; r < informativeReads.size(); r=r+2) {
        RecombReadPair* thisReadPair = new RecombReadPair(informativeReads[r], informativeReads[r+1]);
        informativeReadPairs.push_back(thisReadPair);
    }
    
    
    int numConcordant = 0; int numDiscordant = 0;
    int numMatch = 0; int numMismatch = 0;
    std::vector<double> matchBaseScores; std::vector<double> mismatchBaseScores;
    std::vector<double> concordantBaseScores; std::vector<double> discordantBaseScores;
    std::vector<PhaseSwitch*> phaseSwitches;
    
    int numFullLenghtReadPairs = 0;
    int totalUsedLength = 0;
    std::vector<int> numHets;
    int num0het = 0; int num1het = 0; int num2plusHets = 0;
    
    for (int r = 0; r < informativeReadPairs.size(); r++) {
        
        informativeReadPairs[r]->findAndCombinePairHets(positionToPhase);
        informativeReadPairs[r]->filterHetsByQuality(opt::minBQ);
        
        if (informativeReadPairs[r]->hetSites.size() == 0) {
            num0het++;
           /* std::cout << "goodReadPairs[r]->hetSites.size(): " << goodReadPairs[r]->hetSites.size() << std::endl;
            std::cout << "goodReadPairs[r]->read1->readPos: " << goodReadPairs[r]->read1->readPos << std::endl;
            std::cout << "goodReadPairs[r]->read1->readStrand: " << goodReadPairs[r]->read1->readStrand << std::endl;
            std::cout << "goodReadPairs[r]->read1->usedLength: " << goodReadPairs[r]->read1->usedLength << std::endl;
            std::cout << "goodReadPairs[r]->read1->adjustedReadPos: " << goodReadPairs[r]->read1->adjustedReadPos << std::endl;
            std::cout << "goodReadPairs[r]->read1->CIGAR: " << goodReadPairs[r]->read1->CIGAR << std::endl;
            std::cout << "goodReadPairs[r]->read2->readPos: " << goodReadPairs[r]->read2->readPos << std::endl;
            std::cout << "goodReadPairs[r]->read2->readStrand: " << goodReadPairs[r]->read2->readStrand << std::endl;
            std::cout << "goodReadPairs[r]->read2->usedLength: " << goodReadPairs[r]->read2->usedLength << std::endl;
            std::cout << "goodReadPairs[r]->read2->adjustedReadPos: " << goodReadPairs[r]->read2->adjustedReadPos << std::endl;
            std::cout << "goodReadPairs[r]->read2->CIGAR: " << goodReadPairs[r]->read2->CIGAR << std::endl;
            std::cout <<  std::endl;
            //std::cout << "goodReadPairs[r]->read2->readStrand: " << goodReadPairs[r]->read1-> << std::endl; */
        } else if (informativeReadPairs[r]->hetSites.size() == 1) {
            num1het++;
        } else {
            num2plusHets++;
            if(informativeReadPairs[r]->read1->hetSites.size() >= 1 && informativeReadPairs[r]->read2->hetSites.size() >= 1) {
                *goodReadPairsFile << informativeReadPairs[r]->read1->readPos << "\t" << informativeReadPairs[r]->read2->readPos << std::endl;
            }
        }
        
        /* Just debug
         if (goodReadPairs[r]->read1->CIGAR == "151M" && goodReadPairs[r]->read2->CIGAR == "151M") {
            numFullLenghtReadPairs++;
        } */
        
        totalUsedLength = totalUsedLength + informativeReadPairs[r]->read1->usedLength;
        totalUsedLength = totalUsedLength + informativeReadPairs[r]->read2->usedLength;
        
    }
    
    std::cout << "informativeReadPairs.size(): " << informativeReadPairs.size() << std::endl;
    std::cout << "num0het: " << num0het << std::endl;
    std::cout << "num1het: " << num1het << std::endl;
    std::cout << "num2plusHets: " << num2plusHets << std::endl;
    
    if (opt::hapcutFormat) {
        
        for (int r = 0; r < informativeReadPairs.size(); r++) {
           for (int i = 0; i < informativeReadPairs[r]->hetSites.size(); i++) {
                HetInfo* thisHet = informativeReadPairs[r]->hetSites[i];
                if (thisHet->readPhaseBaseMismatch) {
                    /*  std::cout << "Bases don't match" << std::endl;
                      std::cout << "sequence: " << thisRead->readSeq << std::endl;
                      std::cout << "readBase: " << thisRead->readSeq[i-1] << "[" << readBase << "]" << thisRead->readSeq[i+1] << std::endl;
                      std::cout << "readQual[i]: " << thisRead->readQual[i] << "; score: " << int(thisRead->readQual[i])-33 << std::endl;
                      std::cout << "readPos: " << thisRead->readPos << std::endl;
                      std::cout << "i: " << i << std::endl;
                      std::cout << "readPos + i: " << thisRead->readPos + i << std::endl;
                      std::cout << "phasedSNPbases[0]: " << phasedSNPbases[0] << std::endl;
                      std::cout << "phasedSNPbases[1]: " << phasedSNPbases[1] << std::endl; */
                    numMismatch++;
                    mismatchBaseScores.push_back(thisHet->thisBaseQuality);
                } else {
                    // std::cout << "Fine: " << std::endl;
                     matchBaseScores.push_back(thisHet->thisBaseQuality);
                     numMatch++;
                }
            }
            
            // Now consider all HiC pairs i that map to the chromosome and cover a het at each end - lets call the physical locations of the hets x_i, y_i ( x < y) and let Z_i = 1 if they are in phase, 0 if out of phase (recombined)
            
            // Prod_i ((1-Z_i)(G(y_i)-G(x_i)) + Z_i*(1 -(G(y_i)-G(x_i))
            
            
            
            if (informativeReadPairs[r]->hetSites.size() > 1) {
                int phaseOfPrevious; int posOfPrevious; bool allConcordant = true;
                if (informativeReadPairs[r]->hetSites[0]->thisBase == informativeReadPairs[r]->hetSites[0]->thisHetPhase0) phaseOfPrevious = 0;
                if (informativeReadPairs[r]->hetSites[0]->thisBase == informativeReadPairs[r]->hetSites[0]->thisHetPhase1) phaseOfPrevious = 1;
                for (int i = 1; i < informativeReadPairs[r]->hetSites.size(); i++) {
                    int phaseThis;
                    if (informativeReadPairs[r]->hetSites[i]->thisBase == informativeReadPairs[r]->hetSites[i]->thisHetPhase0) phaseThis = 0;
                    if (informativeReadPairs[r]->hetSites[i]->thisBase == informativeReadPairs[r]->hetSites[i]->thisHetPhase1) phaseThis = 1;
                    
                    if (phaseThis != phaseOfPrevious) {
                        allConcordant = false;
                        PhaseSwitch* thisSwitch = new PhaseSwitch(informativeReadPairs[r]->hetSites[i-1]->pos, informativeReadPairs[r]->hetSites[i]->pos, informativeReadPairs[r]->hetSites[i-1]->thisPhaseQuality, informativeReadPairs[r]->hetSites[i]->thisPhaseQuality);
                        phaseSwitches.push_back(thisSwitch);
                    }
                    posOfPrevious = informativeReadPairs[r]->hetSites[i]->pos;
                    phaseOfPrevious = phaseThis;
                }
                if (allConcordant == false) numDiscordant++;
                else numConcordant++;
            }
        }
        std::cout << "numConcordant: " << numConcordant << std::endl;
        std::cout << "numDiscordant: " << numDiscordant << std::endl;
        
        std::cout << "numMatch: " << numMatch << std::endl;
        std::cout << "numMismatch: " << numMismatch << std::endl;
        std::cout << "Mean mismatchBaseScores: " << vector_average(mismatchBaseScores) << std::endl;
        std::cout << "Mean matchBaseScores: " << vector_average(matchBaseScores) << std::endl;
        
        for (int i = 0; i != phaseSwitches.size(); i++) {
            *phaseSwitchFile << phaseSwitches[i]->posLeft << "\t" << phaseSwitches[i]->posRight << "\t" << phaseSwitches[i]->dist << "\t" << phaseSwitches[i]->phaseQualLeft << "\t" << phaseSwitches[i]->phaseQualRight << std::endl;
        }
        
    } else {
        for (int r = 0; r < informativeReadPairs.size(); r++) {
            if (informativeReadPairs[r]->hetSites.size() > 1) {
                for (int i = 0; i < informativeReadPairs[r]->hetSites.size() - 1; i++) {
                    for (int j = 1; j < informativeReadPairs[r]->hetSites.size(); j++) {
                        char readBase1 = informativeReadPairs[r]->hetSites[i]->thisBase;
                        char readBase2 = informativeReadPairs[r]->hetSites[j]->thisBase;
                        int pos1 = informativeReadPairs[r]->hetSites[i]->pos;
                        int pos2 = informativeReadPairs[r]->hetSites[j]->pos;
                        string pairIndex = numToString(pos1)+"_"+numToString(pos2);
                        if (!SNPpairs.count(pairIndex)) {
                            SNPpairs[pairIndex] = new ReadLinkSNPpair(pos1,pos2,readBase1,readBase2);
                        } else {
                            SNPpairs[pairIndex]->base1.push_back(readBase1);
                            SNPpairs[pairIndex]->base1.push_back(readBase2);
                        }
                    }
                }
            }
        }
        std::cout << "SNPpairs.size(): " << SNPpairs.size() << std::endl;
        std::cout << "numFullLenghtReadPairs: " << numFullLenghtReadPairs << std::endl;
        std::cout << "totalUsedLength: " << totalUsedLength << std::endl;
        std::cout << "totalUsedLength/(goodReadPairs.size()*300): " << (double)totalUsedLength/(double)(informativeReadPairs.size()*300) << std::endl;
    }
    
    return 0;
    
}




void parseDiscordPairsFromSAMOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'b': arg >> opt::minBQ; break;
            case 'm': arg >> opt::minMQ; break;
            case OPT_HAPCUT: opt::hapcutFormat = true; break;
            case 'h':
                std::cout << DISCORDPAIRS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 3) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DISCORDPAIRS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetsFile = argv[optind++];
    opt::samFile = argv[optind++];
}

