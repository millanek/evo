//
//  evo_combineVCFs.h
//  process_vcf
//
//  Created by Milan Malinsky on 27/03/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef evo_combineVCFs_h
#define evo_combineVCFs_h

#include <stdio.h>
#include "process_vcf_utils.h"

void parseVCFcombOptions(int argc, char** argv);
int VCFcombMain(int argc, char** argv);

#endif /* evo_combineVCFs_h */
