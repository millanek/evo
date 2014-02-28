//
//  process_vcf_search_sex.h
//  vcf_process
//
//  Created by Milan Malinsky on 25/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_search_sex__
#define __vcf_process__process_vcf_search_sex__

#include "process_vcf_utils.h"
#include <math.h>

void parseSexSearchOptions(int argc, char** argv);
int sexSearchMain(int argc, char** argv);


class SetDepths {
public:
    SetDepths() : set1_mean_depth(0), set2_mean_depth(0) {};
    
    std::vector<int> set1Depths;
    std::vector<int> set2Depths;
    double set1_mean_depth;
    double set2_mean_depth;
    
    void fillMeanDepths() {
        set1_mean_depth = vector_average(set1Depths);
        set2_mean_depth = vector_average(set2Depths);
    }
};




#endif /* defined(__vcf_process__process_vcf_search_sex__) */
