//
//  evo_diversity_subsampling.h
//  process_vcf
//
//  Created by Milan Malinsky on 13/03/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#ifndef evo_diversity_subsampling_h
#define evo_diversity_subsampling_h

#include "process_vcf_utils.h"
#include "process_vcf_annotation_tools.h"
#include "process_vcf_fst.h"

int subsamplingDxy(int argc, char** argv);
void parseSubsamplingDxyOptions(int argc, char** argv);

#endif /* evo_diversity_subsampling_h */
