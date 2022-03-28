//
//  evo_diversityFromHaps.h
//  process_vcf
//
//  Created by Milan Malinsky on 04/03/2022.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#ifndef evo_diversityFromHaps_h
#define evo_diversityFromHaps_h

#include <stdio.h>
#include "process_vcf_utils.h"
#include "process_vcf_annotation_tools.h"

int regionPiMain(int argc, char** argv);
void parseRegionPiOptions(int argc, char** argv);

#endif /* evo_diversityFromHaps_h */
