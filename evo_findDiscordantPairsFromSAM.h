//
//  evo_findDiscordantPairs.hpp
//  process_vcf
//
//  Created by Milan Malinsky on 08.04.22.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#ifndef evo_findDiscordantPairs_h
#define evo_findDiscordantPairs_h

#include <stdio.h>
#include "evo_recombUtils.h"

void parseDiscordPairsFromSAMOptions(int argc, char** argv);
int DiscordPairsFromSAMMain(int argc, char** argv);


#endif /* evo_findDiscordantPairs_h */
