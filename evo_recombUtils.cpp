//
//  evo_recombUtils.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 29.07.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#include "evo_recombUtils.h"

void RecombRead::findHetsInMatchingString(std::vector<HetInfo*>& hetsOnThisRead, const string& matchSeq, int startPos, const std::map<int,PhaseInfo*>& positionToPhase) {
    for (int i = 0; i < matchSeq.length(); i++) {
        if (positionToPhase.count(readPos + i) == 1) {
            PhaseInfo* thisHetPhase = positionToPhase.at(readPos + i);
            std::vector<char> phasedSNPbases = thisHetPhase->phasedVars;
            char readBase = matchSeq[i];
            int snpPos = readPos + i;
            HetInfo* het = new HetInfo(snpPos, readBase, int(readQual[i])-33, phasedSNPbases[0], phasedSNPbases[1], thisHetPhase->quality);
            hetsOnThisRead.push_back(het);
        }
    }
}

std::vector<HetInfo*> RecombRead::findHetsInRead(const std::map<int,PhaseInfo*>& positionToPhase) {
    std::vector<HetInfo*> hetsOnThisRead;
    
    int startPos = readPos; string processedReadSeq = readSeq;
    while (GIGARtypes.size() > 0) {
        if (GIGARtypes[0] == SOFT_CLIP_CIGAR) {
            readSeq = readSeq.substr(GIGARnums[0]);
            GIGARtypes.erase(GIGARtypes.begin());
            GIGARnums.erase(GIGARnums.begin());
        }

        if (GIGARtypes[0] == HARD_CLIP_CIGAR) {
            GIGARtypes.erase(GIGARtypes.begin());
            GIGARnums.erase(GIGARnums.begin());
        }
        
        if (GIGARtypes[0] == MATCH_CIGAR) {
            string matchSeq = processedReadSeq.substr(0, GIGARnums[0]);
            findHetsInMatchingString(hetsOnThisRead, matchSeq, startPos, positionToPhase);
            startPos = startPos + GIGARnums[0];
            usedLength = usedLength + GIGARnums[0];
            processedReadSeq = processedReadSeq.substr(GIGARnums[0]);
            
            GIGARtypes.erase(GIGARtypes.begin());
            GIGARnums.erase(GIGARnums.begin());
        }
        
        if (GIGARtypes[0] == DELETION_CIGAR) {
            startPos = startPos + GIGARnums[0];
            GIGARtypes.erase(GIGARtypes.begin());
            GIGARnums.erase(GIGARnums.begin());
        }
        
        if (GIGARtypes[0] == INSERTION_CIGAR) {
            processedReadSeq = processedReadSeq.substr(GIGARnums[0]);
            GIGARtypes.erase(GIGARtypes.begin());
            GIGARnums.erase(GIGARnums.begin());
        }
        
    }
        
        
    }
    
    return hetsOnThisRead;
}


string RecombRead::assignStrandFromFlag() {
    string strand = "?";
    if (flag == 81 || flag == 113 || flag == 145 || flag == 177 || flag == 185 || flag == 121) {
        strand = "-";
    } else if (flag == 65 || flag == 73 || flag == 97 || flag == 129 || flag == 161 || flag == 137) {
        strand = "+";
    } else {
        std::cerr << "Unexpected read flag: " << flag << std::endl;
        exit(1);
    }
    return strand;
}

void RecombRead::generateCIGARvectors() {
    string CIGARnum = "";
    for (int i = 0; i < CIGAR.length(); i++) {
       // std::cout << "CIGAR[i]: " << CIGAR[i] << std::endl;
       // std::cout << "isdigit(CIGAR[i]): " << isdigit(CIGAR[i]) << std::endl;
        if (isdigit(CIGAR[i])) {
            CIGARnum += CIGAR[i];
         //   std::cout << "CIGARnum: " << CIGARnum << std::endl;
        } else {
           // std::cout << "CIGARnum: " << CIGARnum << std::endl;
            int CIGARnumInt = atoi(CIGARnum.c_str());
            GIGARnums.push_back(CIGARnumInt);
            if (CIGAR[i] != SOFT_CLIP_CIGAR && CIGAR[i] != INSERTION_CIGAR) GIGARnumsNoSI.push_back(CIGARnumInt);
            GIGARtypes.push_back(CIGAR[i]);
            CIGARnum = "";
        }
    }
}

void RecombReadPair::findAndCombinePairHets(std::map<int,PhaseInfo*>& positionToPhase) {
    read1->hetSites = read1->findHetsInRead(positionToPhase);
    read2->hetSites = read2->findHetsInRead(positionToPhase);
    
    if (read1->hetSites.size() > 0) hetSites = read1->hetSites;
    
    if (read2->hetSites.size() > 0) {
        if (hetSites.size() == 0) hetSites = read2->hetSites;
        else hetSites.insert(hetSites.end(), read2->hetSites.begin(), read2->hetSites.end());
    }
}

void RecombReadPair::filterHetsByQuality(int minQuality) {
    
    std::vector<HetInfo*> goodHets;
    for (std::vector<HetInfo*>::iterator it = hetSites.begin(); it != hetSites.end(); it++) {
        if ((*it)->thisBaseQuality >= minQuality) {
            goodHets.push_back((*it));
        }
    }
    hetSites = goodHets;
}
