//
//  evo_FstAgainstAll.h
//  process_vcf
//
//  Created by Milan Malinsky on 07/02/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef evo_FstAgainstAll_hpp
#define evo_FstAgainstAll_hpp

#include <stdio.h>
#include "process_vcf_utils.h"

void parseFstGlobaloptions(int argc, char** argv);
int FstGlobalMain(int argc, char** argv);

#endif /* evo_FstAgainstAll_h */
