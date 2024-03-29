//
//  evo_recombUtils.h
//  process_vcf
//
//  Created by Milan Malinsky on 29.07.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#ifndef evo_recombUtils_h
#define evo_recombUtils_h

#include <stdio.h>
#include "process_vcf_utils.h"


#define SOFT_CLIP_CIGAR 'S'
#define HARD_CLIP_CIGAR 'H'
#define MATCH_CIGAR 'M'
#define INSERTION_CIGAR 'I'
#define DELETION_CIGAR 'D'

class HetInfo {
    public:
    HetInfo(int position, char readBase, int readBaseQuality, char phase0, char phase1, double phaseQuality, int phaseBlockIn) {
        pos = position;
        thisBase = readBase; thisBaseQuality = readBaseQuality;
        thisHetPhase0 = phase0; thisHetPhase1 = phase1;
        thisPhaseQuality = phaseQuality;
        phaseBlock = phaseBlockIn;
        
        if (thisBase != thisHetPhase0 && thisBase != thisHetPhase1) {
            readPhaseBaseMismatch = true;
        } else {
            readPhaseBaseMismatch = false;
        }
        
        if (thisBase == thisHetPhase0) {
            thisHetPhase01 = 0;
        } else if (thisBase == thisHetPhase1) {
            thisHetPhase01 = 1;
        }
        
    };
    
    int pos;
    int phaseBlock;
    
    char thisBase;
    int thisBaseQuality;
    
    char thisHetPhase0;
    char thisHetPhase1;
    double thisPhaseQuality;
    
    int thisHetPhase01;
    
    bool readPhaseBaseMismatch;

};

class PhaseInfo {
    public:
    PhaseInfo() {};
    
    PhaseInfo(int position, double qual, int cov, std::vector<char>& phasedVarsIn, int blockNumIn) {
        pos = position; quality = qual;
        coverage = cov;
        phasedVars = phasedVarsIn;
        blockNum = blockNumIn;
    };
    
    int pos;
    std::vector<char> phasedVars;
    double quality;
    int coverage;
    int blockNum;
    
};

class PhasePair {
    public:
    PhasePair(): phase1(NULL), phase2(NULL) {};
    
    PhaseInfo* phase1;
    PhaseInfo* phase2;
    
};

class ReadLinkSNPpair {
    public:
    ReadLinkSNPpair(int p1, int p2, char b1, char b2) {
        pos1 = p1; pos2 = p2;
        base1.push_back(b1); base2.push_back(b2);
    };
    
    int pos1; int pos2;
    std::vector<char> base1; std::vector<char> base2;
    
};




class RecombRead {
    public:
    RecombRead(const std::vector<string> samRecVec) : usedLength(0) {
        
        flag = atoi(samRecVec[1].c_str());
        readStrand = assignStrandFromFlag();
        readName = samRecVec[0]; readPos = atoi(samRecVec[3].c_str());
        MQ = atoi(samRecVec[4].c_str()); CIGAR = samRecVec[5];
        readSeq = samRecVec[9]; readQual = samRecVec[10];
        
        generateCIGARvectors();
    };
    
    bool operator< (const RecombRead &other) const {
            return readPos < other.readPos;
    }
    
    
    int flag;
    string readStrand; string readName; int readPos; int adjustedReadPos;
    string readSeq; string readQual;
    int MQ; string CIGAR;
    std::vector<int> GIGARnums; std::vector<char> GIGARtypes;
    std::vector<int> GIGARnumsNoSI;
    
    int usedLength;
    std::vector<HetInfo*> hetSites;
    std::map<int,std::vector<int>> BlockIDsToHetPos;
    
    void findHetsInRead(const std::map<int,PhaseInfo*>& positionToPhase);
    
    private:
        string assignStrandFromFlag();
        void generateCIGARvectors();
        void findHetsInMatchingString(const string& matchSeq, int startPos, const std::map<int,PhaseInfo*>& positionToPhase);

};



class PhaseSwitch {
    public:
    PhaseSwitch(int left, int right, double qLeft, double qRight) {
        posLeft = left;
        posRight = right;
        phaseQualLeft = qLeft;
        phaseQualRight = qRight;
        dist = abs(right - left);
    };
    
    int posLeft;
    int posRight;
    double phaseQualLeft;
    double phaseQualRight;
    int dist;
    

};




class RecombReadPair {
    public:
    RecombReadPair(RecombRead* r1, RecombRead* r2) {
        if (r1 < r2) {
            read1 = r1;
            read2 = r2;
        } else {
            read1 = r2;
            read2 = r1;
        }
    };
    
    bool operator< (const RecombReadPair &other) const {
            return read1->readPos < other.read1->readPos;
    }
    
    
    
    RecombRead* read1;
    RecombRead* read2;
    
    int pairSpan;
    bool pairDiscordant;
    std::vector<HetInfo*> hetSites;
    //std::map<int,std::vector<int>> BlockIDsToHetPos;
    //std::map<int> HetPosToBlockIDs;
    
    void findAndCombinePairHets(const std::map<int,PhaseInfo*> & positionToPhase);
    void filterHetsByQuality(int minQuality);
    void filterHetsByBlock(int blockNum);
    

};

inline unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    
    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

#endif /* evo_recombUtils_h */
