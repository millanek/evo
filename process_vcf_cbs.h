//
//  process_vcf_cbs.h
//  process_vcf
//
//  Created by Milan Malinsky on 22/11/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef __process_vcf__process_vcf_cbs__
#define __process_vcf__process_vcf_cbs__

#include "process_vcf_utils.h"
#include "process_vcf_seq_utils.h"
#include "process_vcf_stats_utils.h"
#include "process_vcf_annotation_tools.h"
#include <unordered_set>

int cbsMain(int argc, char** argv);
void parseCbsOptions(int argc, char** argv);

class cbsSets {
public:
    cbsSets() : set1vsSet2pN(0), sets1and2vsSet3pN(0), withinSet1pN(0), withinSet2pN(0), initialised(false) {};
    
    cbsSets(std::ifstream*& setsFile) {
        string line;
        string set1String; string set2String; string set3String;
        getline(*setsFile, set1String);
        getline(*setsFile, set2String);
        std::vector<string> set1 = split(set1String, ',');
        std::vector<string> set2 = split(set2String, ',');
        for (int i = 0; i < set1.size(); i++) {
            set1Loci.insert(atoi(set1[i].c_str()));
        }
        for (int i = 0; i < set2.size(); i++) {
            set2Loci.insert(atoi(set2[i].c_str()));
        }
        std::cerr << "set 1: "; print_vector(set1, std::cerr);
        std::cerr << "set 2: "; print_vector(set2, std::cerr);
        initialised = true;
    }
    
    std::unordered_set<size_t> set1Loci;
    std::unordered_set<size_t> set2Loci;

    
    double withinSet1pN;
    double withinSet2pN;
    double set1vsSet2pN;
    double sets1and2vsSet3pN;
    bool initialised;
    
};



#endif /* defined(__process_vcf__process_vcf_cbs__) */
