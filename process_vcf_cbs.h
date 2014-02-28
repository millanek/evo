//
//  process_vcf_cbs.h
//  process_vcf
//
//  Created by Milan Malinsky on 22/11/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef __process_vcf__process_vcf_cbs__
#define __process_vcf__process_vcf_cbs__

#include "process_vcf_utils.h"
#include "process_vcf_seq_utils.h"
#include "process_vcf_stats_utils.h"

int cbsMain(int argc, char** argv);
void parseCbsOptions(int argc, char** argv);

#endif /* defined(__process_vcf__process_vcf_cbs__) */
