//
//  process_vcf_mt_sequences.h
//  vcf_process
//
//  Created by Milan Malinsky on 30/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_mt_sequences__
#define __vcf_process__process_vcf_mt_sequences__

#include "process_vcf_utils.h"
#include "process_vcf_seq_utils.h"

void parseGetMtSeqOptions(int argc, char** argv);
int getMtSeqMain(int argc, char** argv);

#endif /* defined(__vcf_process__process_vcf_mt_sequences__) */
