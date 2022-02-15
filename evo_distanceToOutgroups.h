//
//  evo_distanceToOutgroups.hpp
//  process_vcf
//
//  Created by Milan Malinsky on 14/02/2022.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#ifndef evo_distanceToOutgroups_h
#define evo_distanceToOutgroups_h

#include <stdio.h>
#include "process_vcf_utils.h"
#include "process_vcf_annotation_tools.h"


void parseDistOutOptions(int argc, char** argv);
int DistOutMain(int argc, char** argv);

#endif /* evo_distanceToOutgroups_h */
