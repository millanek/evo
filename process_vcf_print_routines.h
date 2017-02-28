//
//  process_vcf_print_routines.h
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef vcf_process_process_vcf_print_routines_h
#define vcf_process_process_vcf_print_routines_h

#include "process_vcf_utils.h"

// Printing doubletons
void print_doubleton_distribution(const string& fileRoot, const std::vector<std::string>& header, std::vector<std::vector<int> >& doubletons);

// Printing het counts
void print_het_counts(const string& fileRoot, const std::vector<std::string>& header, const std::vector<int>& hetCounts, const std::vector<int>& sharedHetCounts);
 

// Printing pairwise difference statistics
void print_pairwise_diff_stats(const string& fileRoot, const std::vector<std::string>& header, const int totalVariantNumber, const std::vector<std::vector<double> >& diffMatrix, const std::vector<std::vector<double> >& diffMatrixMe, const std::vector<std::vector<double> >& diffMatrixHetsVsHomDiff);

void print_H1_pairwise_diff_stats(const string& fileRoot, std::vector<std::string>& header, const int totalVariantNumber, const std::vector<std::vector<double> >& diffMatrixH1);
void print_AllH_pairwise_diff_stats(const string& fileRoot, const std::vector<std::string>& samples, const int totalVariantNumber, const std::vector<std::vector<double> >& diffMatrixAllH);
#endif
