//
//  process_vcf_stats_functions.h
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef vcf_process_process_vcf_stats_functions_h
#define vcf_process_process_vcf_stats_functions_h
#include "process_vcf_utils.h"

// Some helper functions for doing analyses at the level of populations
std::vector<string> initialisePopulationMap(const std::vector<std::string> populations, std::map<std::string,int>& fp_map);
std::map<std::string,std::vector<int>> getPopulationsToIndividualsMap(std::vector<std::string> populations, std::map<std::string,int>& fp_map);

// Increment a counter for each individual who is het at this site
void het_analysis(std::vector<int>& hetCounts, std::vector<int>& sharedHetCounts, FilterResult& result);


// DOUBLETONS:
// Initializing the doubleton data structure if counting doubletons per population
void initializeDoubletons(std::vector<std::vector<int> >& d, const std::vector<std::string> populations);
// Performing doubleton analysis
void doubleton_analysis(std::vector<std::vector<int> >& doubletons, FilterResult& result, int numChromosomes, const std::vector<std::string> p, std::map<std::string,int>& p_int_map);
void privateVars_analysis(std::vector<int>& privateVarCounts, const FilterResult& result, const std::vector<std::vector<size_t> >& populationsIndices, const std::vector<std::vector<size_t> >& populationsIndicesComplements);

// Increment the diff (half)matrix, recording the differences 
void diffs_between_individuals(std::vector<std::vector<double> >& diffs,std::vector<std::vector<double> >& diffs_me, std::vector<std::vector<double> >& diffs_Hets_vs_Homs, std::vector<std::vector<int> >& pairwise_missingness, FilterResult& result);

// Increment the diff (half)matrix, using haplotypes
void diffs_between_H1(std::vector<std::vector<double> >& H1diffs, FilterResult& result);
void diffs_between_AllH(std::vector<std::vector<double> >& AllHdiffs, FilterResult& result);
#endif
