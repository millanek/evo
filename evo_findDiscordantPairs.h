//
//  evo_findDiscordantPairs.hpp
//  process_vcf
//
//  Created by Milan Malinsky on 08.04.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#ifndef evo_findDiscordantPairs_h
#define evo_findDiscordantPairs_h

#include <stdio.h>
#include "process_vcf_utils.h"

#define SOFT_CLIP_CIGAR 'S'
#define HARD_CLIP_CIGAR 'H'
#define MATCH_CIGAR 'M'
#define INSERTION_CIGAR 'I'

void parseDiscordPairsOptions(int argc, char** argv);
int DiscordPairsMain(int argc, char** argv);


class HetInfo {
    public:
    HetInfo(int position, char readBase, int readBaseQuality, char phase0, char phase1, double phaseQuality) {
        pos = position;
        thisBase = readBase; thisBaseQuality = readBaseQuality;
        thisHetPhase0 = phase0; thisHetPhase1 = phase1;
        thisPhaseQuality = phaseQuality;
        
        if (thisBase != thisHetPhase0 && thisBase != thisHetPhase1) {
            readPhaseBaseMismatch = true;
        } else {
            readPhaseBaseMismatch = false;
        }
    };
    
    int pos;
    
    char thisBase;
    int thisBaseQuality;
    
    char thisHetPhase0;
    char thisHetPhase1;
    double thisPhaseQuality;
    
    bool readPhaseBaseMismatch;

};

class PhaseInfo {
    public:
    
    PhaseInfo(int position, double qual, int cov, std::vector<char>& phasedVarsIn) {
        pos = position; quality = qual;
        coverage = cov;
        phasedVars = phasedVarsIn;
    };
    
    int pos;
    std::vector<char> phasedVars;
    double quality;
    int coverage;
    
};

class RecombRead {
    public:
    RecombRead(const std::vector<string> samRecVec) {
        
        flag = atoi(samRecVec[1].c_str());
        readStrand = assignStrandFromFlag();
        readName = samRecVec[0]; readPos = atoi(samRecVec[3].c_str());
        MQ = atoi(samRecVec[4].c_str()); CIGAR = samRecVec[5];
        readSeq = samRecVec[9]; readQual = samRecVec[10];
        
        generateCIGARvectors();
    };
    
    int flag;
    string readStrand; string readName; int readPos;
    string readSeq; string readQual;
    int MQ; string CIGAR;
    std::vector<int> GIGARnums; std::vector<char> GIGARtypes;
    std::vector<int> GIGARnumsNoSI;
    
    std::vector<HetInfo*> findHetsInRead(std::map<int,PhaseInfo*>& positionToPhase);
    
    private:
        string assignStrandFromFlag();
        void generateCIGARvectors();

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
        read1 = r1;
        read2 = r2;
    };
    
    RecombRead* read1;
    RecombRead* read2;
    
    int pairSpan;
    bool pairDiscordant;
    std::vector<HetInfo*> hetSites;
    
    void findAndCombinePairHets(std::map<int,PhaseInfo*>& positionToPhase);
    void filterHetsByQuality(int minQuality);
    

};

#endif /* evo_findDiscordantPairs_h */
