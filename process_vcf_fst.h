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


inline double calculateFstNumerator(const double p1, const double p2, const int n1, const int n2) {
    double power = pow((p1-p2), 2);
    double fraction1 = (p1*(1-p1))/(n1-1);
    double fraction2 = (p2*(1-p2))/(n2-1);
    double numerator = power - fraction1 - fraction2;
    return numerator;
}

inline double calculateFstDenominator(const double p1, const double p2) {
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

// For now this requires that the VCF has one individual per species
inline double calculateOverallDxy(const Counts& thisVarCounts) {
    double Dxy;
    int numSamples = (int)thisVarCounts.individualsWithVariant.size();
    int sumKij = 0;
    for (std::vector<std::string>::size_type i = 0; i != numSamples - 1; i++) {
        for (std::vector<std::string>::size_type j = i+1; j != numSamples; j++) {
            if (thisVarCounts.individualsWithVariant[i] == 0 && thisVarCounts.individualsWithVariant[j] == 0) {
            } else if (thisVarCounts.individualsWithVariant[i] == 1 && thisVarCounts.individualsWithVariant[j] == 0) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.individualsWithVariant[i] == 0 && thisVarCounts.individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.individualsWithVariant[i] == 1 && thisVarCounts.individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.individualsWithVariant[i] == 2 && thisVarCounts.individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.individualsWithVariant[i] == 1 && thisVarCounts.individualsWithVariant[j] == 2) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.individualsWithVariant[i] == 2 && thisVarCounts.individualsWithVariant[j] == 2) {
            } else if (thisVarCounts.individualsWithVariant[i] == 2 && thisVarCounts.individualsWithVariant[j] == 0) {
                sumKij = sumKij + 4;
            } else if (thisVarCounts.individualsWithVariant[i] == 0 && thisVarCounts.individualsWithVariant[j] == 2) {
                sumKij = sumKij + 4;
            }
        }
    }
    Dxy = (double)sumKij/(2*(numSamples*(numSamples-1)));
    return Dxy;
}



#endif /* defined(__process_vcf__process_vcf_fst__) */
