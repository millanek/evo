//
//  process_vcf_polymorphic_in_sets.h
//  process_vcf
//
//  Created by Milan Malinsky on 31/10/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef __process_vcf__process_vcf_polymorphic_in_sets__
#define __process_vcf__process_vcf_polymorphic_in_sets__

#include "process_vcf_utils.h"

int polymorphicMain(int argc, char** argv);
void parsePolymorphicOptions(int argc, char** argv);


std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);

#endif /* defined(__process_vcf__process_vcf_polymorphic_in_sets__) */
