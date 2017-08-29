//
//  evo_shared_variation.h
//  vcf_process
//
//  Created by Milan Malinsky on 25/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_fixed_search__
#define __vcf_process__process_vcf_fixed_search__

#include "process_vcf_utils.h"

int sharedVarMain(int argc, char** argv);
void parseSharedVarOptions(int argc, char** argv);


std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);


#endif /* defined(__vcf_process__process_vcf_fixed_search__) */
