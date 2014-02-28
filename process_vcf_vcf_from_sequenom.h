//
//  process_vcf_vcf_from_sequenom.h
//  process_vcf
//
//  Created by Milan Malinsky on 25/10/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef __process_vcf__process_vcf_vcf_from_sequenom__
#define __process_vcf__process_vcf_vcf_from_sequenom__

#include "process_vcf_utils.h"

void parseVCFfromSequenomOptions(int argc, char** argv);
int VCFfromSequenomMain(int argc, char** argv);

static const char *VCF_HEADER =
"##fileformat=VCFv4.1\n"
"##source=SequenomGenotyping\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

#endif /* defined(__process_vcf__process_vcf_vcf_from_sequenom__) */
