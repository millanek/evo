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
    
    std::map<int,PhaseInfo*> posToPhase;
    std::map<string,std::vector<int>> infoPairNameToPos;
    std::map<string,std::vector<string>> infoPairNameToStrands;
    std::vector<int> phaseBlockSNPnums;
    
    std::map<string,std::vector<RecombRead*>> samNameToReads;
    int numHetPairs = 0;
    std::unordered_map<string, ReadLinkSNPpair*> SNPpairs;
    
    if (opt::hapcutFormat) {
        int blockNum = 0;
        // Parse the Hapcut blocks file
        while (getline(*hetsFile, line)) {
            if (line[0] == '*') {
            
            } else if (line[0] == 'B' && line[1] == 'L') { // New block - should in the future separate the hets by blocks
                blockNum++; phaseBlockSNPnums.push_back(0);
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
                PhaseInfo* thisPhase = new PhaseInfo(snpPos,phaseQual,snpCoverage, phasedVars,blockNum);
                posToPhase[snpPos] = thisPhase;
                phaseBlockSNPnums[blockNum-1]++;
            }
        }
        //maxBlockIndex = (int)std::distance(phaseBlockSNPnums.begin(),std::max_element(phaseBlockSNPnums.begin(), phaseBlockSNPnums.end()));
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
            PhaseInfo* thisPhase = new PhaseInfo(snpPos,phaseQual,snpCoverage, phasedVars, 1);
            posToPhase[snpPos] = thisPhase;
        }
        numHetPairs = nChoosek((int)posToPhase.size(),2);
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
    
    int num0het = 0; int num1het = 0; int num2plusHets = 0;
    int totalUsedLength = 0;
    std::vector<RecombReadPair*> informativeReadPairs;
    for (int r = 0; r < informativeReads.size(); r=r+2) {
        RecombReadPair* thisReadPair = new RecombReadPair(informativeReads[r], informativeReads[r+1]);
        thisReadPair->findAndCombinePairHets(posToPhase);
        thisReadPair->filterHetsByQuality(opt::minBQ);
        
        if (thisReadPair->hetSites.size() == 0) {
            num0het++;
        } else if (thisReadPair->hetSites.size() == 1) {
            num1het++;
        } else {
            num2plusHets++;
        }
        
        if (thisReadPair->hetSites.size() > 1) {
            for (std::map<int, std::vector<int>>::iterator it = thisReadPair->read1->BlockIDsToHetPos.begin();
                 it != thisReadPair->read1->BlockIDsToHetPos.end(); it++) {
                if (thisReadPair->read2->BlockIDsToHetPos.count(it->first) == 1) {
                    informativeReadPairs.push_back(thisReadPair);
                    totalUsedLength += thisReadPair->read1->usedLength;
                    totalUsedLength += thisReadPair->read2->usedLength;
                }
            }
        }
    }
    
    std::cout << "Initial Read Pairs.size(): " << informativeReads.size()/2.0 << std::endl;
    std::cout << "num0het: " << num0het << std::endl;
    std::cout << "num1het: " << num1het << std::endl;
    std::cout << "num2plusHets: " << num2plusHets << std::endl;
    std::cout << "informativeReadPairs.size(): " << informativeReadPairs.size() << std::endl;
    
    
    
    int numConcordant = 0; int numDiscordant = 0;
    int numMatch = 0; int numMismatch = 0; int totalEffectiveLength = 0;
    std::vector<double> matchBaseScores; std::vector<double> mismatchBaseScores;
    std::vector<double> concordantBaseScores; std::vector<double> discordantBaseScores;
    std::vector<PhaseSwitch*> phaseSwitches;
    std::vector<std::vector<int>> phaseConcordanceCoords;
        
    
    if (opt::hapcutFormat) {
        for (int r = 0; r < informativeReadPairs.size(); r++) {
            
           for (int i = 0; i < informativeReadPairs[r]->hetSites.size(); i++) {
                HetInfo* thisHet = informativeReadPairs[r]->hetSites[i];
                if (thisHet->readPhaseBaseMismatch) {
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
            
            // std::vector<PhaseSwitch*> thisPairSwitches;
            std::vector<int> switchPairI; std::vector<int> switchPairJ;
            std::vector<int> concordPairI; std::vector<int> concordPairJ;
            for (int i = 0; i < informativeReadPairs[r]->hetSites.size() - 1; i++) {
                for (int j = 1; j < informativeReadPairs[r]->hetSites.size(); j++) {
                    if (informativeReadPairs[r]->hetSites[i]->phaseBlock == informativeReadPairs[r]->hetSites[j]->phaseBlock) {
                        int phaseI = informativeReadPairs[r]->hetSites[i]->thisHetPhase01;
                        int phaseJ = informativeReadPairs[r]->hetSites[j]->thisHetPhase01;
                        if (phaseI != phaseJ) {
                            switchPairI.push_back(i); switchPairJ.push_back(j);
                            //int iPos = informativeReadPairs[r]->hetSites[i]->pos;
                            //int jPos = informativeReadPairs[r]->hetSites[j]->pos;
                            //int iQual = informativeReadPairs[r]->hetSites[i]->thisPhaseQuality;
                            //int jQual = informativeReadPairs[r]->hetSites[j]->thisPhaseQuality;
                            // PhaseSwitch* thisSwitch = new PhaseSwitch(iPos, jPos, iQual, jQual);
                            // thisPairSwitches.push_back(thisSwitch);
                        } else {
                            concordPairI.push_back(i); concordPairJ.push_back(j);
                        }
                    }
                }
            }
            
            // TO DO:
            // Select the 'right' switch pair if there are multiple options:
            // The shortest one? Needs more thought....
            if (switchPairI.size() > 0) {
                numDiscordant++;
                int iPos = informativeReadPairs[r]->hetSites[switchPairI[0]]->pos;
                int jPos = informativeReadPairs[r]->hetSites[switchPairJ[0]]->pos;
                int iQual = informativeReadPairs[r]->hetSites[switchPairI[0]]->thisPhaseQuality;
                int jQual = informativeReadPairs[r]->hetSites[switchPairJ[0]]->thisPhaseQuality;
                PhaseSwitch* thisSwitch = new PhaseSwitch(iPos, jPos, iQual, jQual);
                phaseSwitches.push_back(thisSwitch);
                switchPairI.empty(); switchPairJ.empty();
                totalEffectiveLength = totalEffectiveLength + (jPos - iPos);
            } else {
                numConcordant++;
                std::vector<int> thisConcordantCoords;
                int maxD = 0; int maxDindex = 0;
                for (int i = 0; i != concordPairI.size(); i++) {
                    int iPos = informativeReadPairs[r]->hetSites[concordPairI[i]]->pos;
                    int jPos = informativeReadPairs[r]->hetSites[concordPairJ[i]]->pos;
                    if (jPos - iPos > maxD) {
                        maxDindex = i;
                    }
                }
                int iPosDindex = informativeReadPairs[r]->hetSites[concordPairI[maxDindex]]->pos;
                int jPosDindex = informativeReadPairs[r]->hetSites[concordPairJ[maxDindex]]->pos;
                totalEffectiveLength = totalEffectiveLength + (jPosDindex - iPosDindex);
                std::cout << "read r: " << r << std::endl;
                std::cout << "iPosDindex: " << iPosDindex << std::endl;
                std::cout << "jPosDindex: " << jPosDindex << std::endl;
                
                thisConcordantCoords.push_back(iPosDindex);
                thisConcordantCoords.push_back(jPosDindex);
                phaseConcordanceCoords.push_back(thisConcordantCoords);
            }
            /*if (readPairsProcessed % 100 == 0) {
                std::cout << "readPairsProcessed: " << readPairsProcessed << std::endl;
                std::cout << "informativeReadPairs[r]->hetSites.size(): " << informativeReadPairs[r]->hetSites.size() << std::endl;
                std::cout << "maxHetsNum: " << maxHetsNum << std::endl;
                std::cout << "phaseSwitches.size(): " << phaseSwitches.size() << std::endl;
            } */
        }
        
        std::cout << "Effective coverage (bp): " << totalEffectiveLength << std::endl;
        std::cout << "numConcordant: " << numConcordant << std::endl;
        std::cout << "numDiscordant: " << numDiscordant << std::endl;
        std::cout << "phaseConcordanceCoords.size(): " << phaseConcordanceCoords.size() << std::endl;
        
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
        std::cout << "\n" << DISCORDPAIRS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetsFile = argv[optind++];
    opt::samFile = argv[optind++];
}

