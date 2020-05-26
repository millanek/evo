//
//  evo_AlleleFeq.h
//  process_vcf
//
//  Created by Milan Malinsky on 26/05/2020.
//  Copyright © 2020 Milan Malinsky. All rights reserved.
//

#ifndef evo_AlleleFeq_h
#define evo_AlleleFeq_h

#include <stdio.h>
#include "process_vcf_utils.h"

void parseAFoptions(int argc, char** argv);
int AFmain(int argc, char** argv);

#endif /* evo_AlleleFeq_h */
