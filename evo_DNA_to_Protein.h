//
//  evo_DNA_to_Protein.h
//  process_vcf
//
//  Created by Milan Malinsky on 01/02/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#ifndef evo_DNA_to_Protein_h
#define evo_DNA_to_Protein_h

#include "process_vcf_utils.h"
#include "process_vcf_coding_sequences.h"

void parseDNAtoProteinOptions(int argc, char** argv);
int DNAtoProtein(int argc, char** argv);

#endif /* evo_DNA_to_Protein_h */
