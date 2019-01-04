//
//  evo_PBS.hpp
//  process_vcf
//
//  Created by Milan Malinsky on 03/01/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef evo_PBS_h
#define evo_PBS_h

#include <stdio.h>
#include "process_vcf_utils.h"

void parsePBSoptions(int argc, char** argv);
int PBSmain(int argc, char** argv);

#endif /* evo_PBS_h */
