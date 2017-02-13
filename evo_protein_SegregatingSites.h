//
//  evo_protein_SegregatingSites.h
//  process_vcf
//
//  Created by Milan Malinsky on 03/02/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#ifndef evo_protein_SegregatingSites_h
#define evo_protein_SegregatingSites_h

#include "process_vcf_utils.h"
#include "process_vcf_coding_sequences.h"

void parseProteinSsOptions(int argc, char** argv);
int ProteinSs(int argc, char** argv);

#endif /* evo_protein_SegregatingSites_h */
