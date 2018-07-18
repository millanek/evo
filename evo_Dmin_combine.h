//
//  evo_Dmin_combine.h
//  process_vcf
//
//  Created by Milan Malinsky on 13/07/2018.
//  Copyright © 2018 Milan Malinsky. All rights reserved.
//

#ifndef evo_Dmin_combine_h
#define evo_Dmin_combine_h

#include "process_vcf_utils.h"
#include "process_vcf_stats_utils.h"

void parseDminCombineOptions(int argc, char** argv);
int DminCombineMain(int argc, char** argv);

#endif /* evo_Dmin_combine_h */
