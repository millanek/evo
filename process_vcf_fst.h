//
//  process_vcf_fst.h
//  process_vcf
//
//  Created by Milan Malinsky on 03/11/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef __process_vcf__process_vcf_fst__
#define __process_vcf__process_vcf_fst__

#include "process_vcf_utils.h"
#include <math.h>

void parseFstOptions(int argc, char** argv);
int fstMain(int argc, char** argv);


inline double calculateFstNumerator(const SetCounts& thisVarCounts, const int n1, const int n2) {
    double p1 = (double)thisVarCounts.set1Count/n1;
    double p2 = (double)thisVarCounts.set2Count/n2;
    double power = pow((p1-p2), 2);
    double fraction1 = (p1*(1-p1))/(n1-1);
    double fraction2 = (p2*(1-p2))/(n2-1);
    double numerator = power - fraction1 - fraction2;
    return numerator;
}

inline double calculateFstDenominator(const SetCounts& thisVarCounts, const int n1, const int n2) {
    double p1 = (double)thisVarCounts.set1Count/n1;
    double p2 = (double)thisVarCounts.set2Count/n2;
    double denominator = (p1*(1-p2))+(p2*(1-p1));
    return denominator;
}


inline double calculateExpectedHeterozygositySimple(const double p) {
    double q = 1 - p;
    double heterozygosity = 1 - (pow(p,2)+pow(q,2));
    return heterozygosity;
}


inline double calculateExpectedHeterozygosityNei78(const double p, const int n) {
    double q = 1 - p;
    double simpleHeterozygosity = 1 - (pow(p,2)+pow(q,2));
    double heterozygosity = (n*simpleHeterozygosity)/(n-1);
    return heterozygosity;
}

#endif /* defined(__process_vcf__process_vcf_fst__) */
