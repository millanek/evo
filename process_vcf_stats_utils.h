//
//  process_vcf_stats_utils.h
//  vcf_process
//
//  Created by Milan Malinsky on 01/10/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_stats_utils__
#define __vcf_process__process_vcf_stats_utils__

#include "process_vcf_utils.h"
#include <math.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace boost::math;
using boost::math::chi_squared;
using boost::math::students_t;
using boost::math::quantile;
using boost::math::complement;

// -------------------------------------    BASIC MATH/STATS  ----------------------------------------

// factorial(x): (x! for non-negative integer x) is defined to be gamma(x+1) (as in R)
inline double factorial(double num) {
    if (num < 0) {
        std::cerr << "Can't compute factorial of a negative number " << num << std::endl;
        exit(1);
    }
    return tgamma(num+1);
}

// Calculates the binomial coefficient (n choose k)
inline int choose(int n, int k) {
    double dResult = factorial(n)/(factorial(k)*factorial(n-k));
    int iResult = (int)round(dResult);
    return iResult;
}

// standard deviation of a sample
template <class T> double std_dev(T vector) {
    double mean = vector_average(vector);
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += pow((vector[i] - mean), 2.0);
    }
    double std_dev = sqrt((double)sum/(double)(vector.size()-1));
    return std_dev;
}


inline void copy_except(int i, std::vector<double>& inVec, std::vector<double>& outVec) {
    std::copy(inVec.begin(), inVec.begin() + i, outVec.begin());
    std::copy(inVec.begin() + i + 1, inVec.end(), outVec.begin()+i);
    //std::cerr << "copying:" << i << " "; print_vector_stream(inVec, std::cerr);
    //std::cerr << "copied: " << i << " "; print_vector_stream(outVec, std::cerr);
}

// jackknive standard error
template <class T> double jackknive_std_err(T& vector) {
    std::vector<double> jackkniveAverages;
    std::vector<double> JregionDs; JregionDs.resize(vector.size()-1);
    for (std::vector<double>::size_type i = 0; i != vector.size(); i++) {
        // std::cerr << "copying " << i << std::endl;
        copy_except(i, vector, JregionDs);
        jackkniveAverages.push_back(vector_average(JregionDs));
        JregionDs.clear(); JregionDs.resize(vector.size()-1);
    }
    double jackkniveOverallMean = vector_average(jackkniveAverages);
    double sum = 0;
    for (int i = 0; i < jackkniveAverages.size(); i++) {
        sum += pow((jackkniveAverages[i] - jackkniveOverallMean), 2.0);
    }
    double var = ((double)(jackkniveAverages.size()-1)/(double)jackkniveAverages.size()) * sum;
    double Dstd_err = sqrt(var);
    return Dstd_err;
}

// ---------------------------------------    DENSITY FUNCTIONS  -------------------------------------

// Central chi-squared probability density function (as in R)
// f_n(x) = 1 / (2^(n/2) Γ(n/2)) x^(n/2-1) e^(-x/2)
inline double chisq_pdf(double df, double chi_sq) {
    double denominator = pow(2.0, (df/2.0)) * tgamma(df/2.0) * pow(chi_sq, (df/2.0)-1) * exp(-chi_sq/2.0);
    double density = 1.0/denominator;
    return density;
}

// Get chi-sq pvalues
inline double chisq_cdf(double df, double chi_sq) {
    chi_squared dist(df);
    double p = cdf(dist, chi_sq);
    return p;
}

// Get student's t p-values
inline double students_t_cdf(double df, double t) {
    students_t dist(df);
    double p = cdf(dist, t);
    return p;
}


// ---------------------------------------------    TESTS   ------------------------------------------

// Welch's Two Sample t-test
template <class T> double two_sample_t(T vector1, T vector2, double d = 0) {
    double variance1 = pow(std_dev(vector1), 2.0);
    double variance2 = pow(std_dev(vector2), 2.0);
    double SE = sqrt((variance1/vector1.size())+(variance2/vector2.size()));
    double t = ((vector_average(vector1)-vector_average(vector2)) - d)/SE;
    // DF = (s1^2/n1 + s2^2/n2)^2/{ [ (s1^2 / n1)^2 / (n1 - 1) ] + [ (s2^2 / n2)^2 / (n2 - 1) ] }
    double DFnumerator = pow((variance1/vector1.size())+(variance2/vector2.size()), 2.0);
    double DFdenominatorPt1 = pow(variance1/vector1.size(), 2.0) / (vector1.size() - 1);
    double DFdenominatorPt2 = pow(variance2/vector2.size(), 2.0) / (vector2.size() - 1);
    double DF = DFnumerator/(DFdenominatorPt1 + DFdenominatorPt2);
    
   // std::cerr << "t-stat: " << t << std::endl;
   // std::cerr << "DF: " << DF << std::endl;
    
    if (d == 0) {
        // alternative hypothesis: true difference in means is not equal to 0
        return 2*students_t_cdf(DF, t);
    } else {
        // alternative hypothesis: true difference in means is greater than d
        return 1-students_t_cdf(DF, t);
    }
    
}


// Pearson's chi-squared test of independence
// (a - row1,column1; b - row1,column2; c - row2,column1; d- row2,column2)
inline double pearson_chi_sq_indep(int a, int b, int c, int d) {
    // Expected values for a,b,c,and d: (RowTotal*ColTotal)/GridTotal
    double n = (double)a+b+c+d;
    double exp_a = (double)((a+b)*(a+c))/n;
    double exp_b = (double)((a+b)*(b+d))/n;
    double exp_c = (double)((c+d)*(a+c))/n;
    double exp_d = (double)((c+d)*(b+d))/n;
    //std::cerr << "exp_a: " << exp_a << "exp_b: " << exp_b << "exp_c: " << exp_c << "exp_d: " << exp_d << std::endl;
    double chi_sq = pow((a-exp_a),2.0)/exp_a + pow((b-exp_b),2.0)/exp_b + pow((c-exp_c),2.0)/exp_c + pow((d-exp_d),2.0)/exp_d;
    //std::cerr << "Chisq stat: " << chi_sq << std::endl;
    double df = 1.0;
    double p = (1-chisq_cdf(df, chi_sq));
    return p;
}


// Fisher's exact test: (a - row1,column1; b - row1,column2; c - row2,column1; d- row2,column2)
// Probability of a specific table
inline double fisher_exactTable(int a, int b, int c, int d) {
    
    long double numerator = factorial(a+b) * factorial(c+d) * factorial(a+c) * factorial(b+d);
    long double denominator = factorial(a) * factorial(b) * factorial(c) * factorial(d) * factorial(a+b+c+d);
    
    double p = numerator/denominator;
    return p;
}

// Fisher's exact test: (a - row1,column1; b - row1,column2; c - row2,column1; d- row2,column2)
// Two-tailed p-value
inline double fisher_exact(int a, int b, int c, int d) {
    
    std::vector<int> mTotals;
    int r1 = a + b; int r2 = c + d;
    int c1 = a + c; int c2 = b + d;
    mTotals.push_back(r1); mTotals.push_back(r2); mTotals.push_back(c1); mTotals.push_back(c2);
    int m = std::min(std::min(r1, r2), std::min(c1,c2)); // The smallest margin total
    // int min_index = (int)(std::min_element(mTotals.begin(), mTotals.end()) - mTotals.begin()); // Position of the ssmallest margin total
    
    std::vector<double> allPvals;
    if (m == c2) {
        for (int i = 0; i <= m; i++) {
            double thisP = fisher_exactTable(r1-i,i,r2-(m-i),m-i);
            allPvals.push_back(thisP);
        }
    } else if (m == c1) {
        for (int i = 0; i <= m; i++) {
            double thisP = fisher_exactTable(i,r1-i,m-i,r2-(m-i));
            allPvals.push_back(thisP);
        }
    } else if (m == r2) {
        for (int i = 0; i <= m; i++) {
            double thisP = fisher_exactTable(c1-(m-i),c2-i,m-i,i);
            allPvals.push_back(thisP);
        }
    } else if (m == r1){
        for (int i = 0; i <= m; i++) {
            double thisP = fisher_exactTable(m-i,i,c1-(m-i),c2-i);
            allPvals.push_back(thisP);
        }
    } else {
        std::cerr << "a: " << a << " b: " << b << std::endl;
        std::cerr << "c: " << c << " d: " << d << std::endl;
        std::cerr << "m: " << m << " c1: " << c1 << " c2: " << c2 << " r1: " << r1 << " r2: " << r2 << std::endl;
        assert(false);
    }
    
    double thisTablePval = fisher_exactTable(a,b,c,d);
    double p = 0;
    for (std::vector<double>::size_type i = 0; i != allPvals.size(); i++) {
        if (allPvals[i] <= thisTablePval) {
            p += allPvals[i];
        }
    }
    /*if (m == r2) {
        std::cerr << "a: " << a << " b: " << b << std::endl;
        std::cerr << "c: " << c << " d: " << d << std::endl;
        std::cerr << "p: " << p << std::endl;
    } */
    
//    std::cerr << "this pval:" << thisTablePval << std::endl;
//    print_vector_stream(allPvals, std::cerr);
    return p;
}





#endif /* defined(__vcf_process__process_vcf_stats_utils__) */
