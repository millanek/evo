//
//  evo_getInformativeReadsFromSam.h
//  process_vcf
//
//  Created by Milan Malinsky on 05.04.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#ifndef evo_getInformativeReadsFromSam_h
#define evo_getInformativeReadsFromSam_h

#include <stdio.h>
#include "process_vcf_utils.h"
#include "evo_recombUtils.h"

void parseInfoReadsOptions(int argc, char** argv);
int InfoReadsMain(int argc, char** argv);

#endif /* evo_getInformativeReadsFromSam_h */
