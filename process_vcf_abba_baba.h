//
//  process_vcf_abba_baba.h
//  process_vcf
//
//  Created by Milan Malinsky on 12/12/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef __process_vcf__process_vcf_abba_baba__
#define __process_vcf__process_vcf_abba_baba__

#include "process_vcf_utils.h"

class ABBA_BABA_Freq_allResults {
public:
    ABBA_BABA_Freq_allResults() : Dnumerator(0), Ddenominator(0), lastVarsDnum(0), lastVarsDdenom(0), windowDnum(0), windowDdenom(0), f_d_denominator(0), window_f_d_denominator(0),  lastVarsF_d_denom(0), f_G_denom(0), f_G_num(0), lastVarsF_G_num(0), lastVarsF_G_denom(0) {};
    
    double Dnumerator; double Ddenominator;     // simple D statistic
    double lastVarsDnum; double lastVarsDdenom; // D within a long stretch window for jackkinive analysis
    double windowDnum; double windowDdenom; // D within a window
    double f_d_denominator; double window_f_d_denominator; double lastVarsF_d_denom;
    double f_G_denom; double f_G_num; double lastVarsF_G_num; double lastVarsF_G_denom;
};

int abbaBabaMain(int argc, char** argv);
void parseAbbaBabaOptions(int argc, char** argv);

#endif /* defined(__process_vcf__process_vcf_abba_baba__) */
