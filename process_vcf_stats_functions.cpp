//
//  process_vcf_stats_functions.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_stats_functions.h"

static const double DIFF_WEIGHT_BOTH_HETS_ME = 0;
static const double DIFF_WEIGHT_BOTH_HETS_RICHARD = (2.0/3.0);
static const double DIFF_WEIGHT_HOM_DIFFERENCE_ME = 1;
static const double DIFF_WEIGHT_HOM_DIFFERENCE_RICHARD = (2.0/3.0);
static const double DIFF_WEIGHT_ONE_HOM_ONE_HET = 0.5;

// Increment a counter for each individual who is het at this site
void het_analysis(std::vector<int>& hetCounts, FilterResult& result) {
    for (std::vector<int>::size_type i = 0; i < hetCounts.size(); i++) {
        if (result.counts.individualsWithVariant[i] == 1) {
            hetCounts[i]++;
        }
    }
}


// Initializing the doubleton data structure
std::vector<std::string> initializeDoubletons(std::vector<std::vector<int> >& d, const std::vector<std::string> populations, std::map<std::string,int>& fp_map) {
    std::vector<std::string> pop_unique = populations;
    //std::sort(pop_unique.begin(), pop_unique.end());
    //std::vector<std::string>::iterator it = std::unique (pop_unique.begin(), pop_unique.end());
    //pop_unique.resize(std::distance(pop_unique.begin(), it));
    std::vector<std::string>::size_type n_pop = pop_unique.size();
    
    // Initialize doubletons
    for (int i = 0; i < n_pop; i++) { 
        std::vector<int> v(n_pop,0);
        d.push_back(v);
    }   
    // Initialize fields to populations map
    for (int i = 0; i < n_pop; i++) { 
        fp_map[pop_unique[i]] = i;
    }
//#ifdef TESTING
    std::cerr << "Doubletons initialised with size: " << n_pop << std::endl;
//#endif
    return pop_unique;
}

// Initializing the doubleton data structure to get doubleton distribution per individual
void initializeDoubletonsIndividuals(std::vector<std::vector<int> >& d, const std::vector<std::string> individuals) {
    std::vector<std::string>::size_type n = individuals.size();
    
    // Initialize doubletons
    for (int i = 0; i < n; i++) { 
        std::vector<int> v(n,0);
        d.push_back(v);
    }   
}


// Filling in the doubleton data structure
void doubleton_analysis(std::vector<std::vector<int> >& doubletons, FilterResult& result, int numChromosomes, const std::vector<std::string> p, std::map<std::string,int>& p_int_map) {
//#ifdef TESTING
    std::cerr << "Overall count: " << result.counts.overall << std::endl;
//#endif
    int fields[2];
    
    if (result.counts.overall == 2) {
        std::vector<int>::iterator it_hom;
        std::vector<int>::iterator it_het;
        // Check if the doubleton is in a single homozygous individual
        it_hom = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 2);
        if (it_hom != result.counts.individualsWithVariant.end()) {
            int doubleton_hom = std::distance(result.counts.individualsWithVariant.begin(), it_hom);
            fields[0] = doubleton_hom; fields[1] = doubleton_hom;
        } else {
            // Otherwise look at which two individuals have the two heterozygous variants
            it_het = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 1);
            int i = 0;
            while (it_het != result.counts.individualsWithVariant.end()) {
                int doubleton_het = std::distance(result.counts.individualsWithVariant.begin(), it_het);
                fields[i] = doubleton_het;
                ++it_het;
                it_het = std::find(it_het, result.counts.individualsWithVariant.end(), 1);
                ++i;
            }
            assert(i <= 3);
        }
#ifdef TESTING
        std::cerr << "Fields 0: " << fields[0] << " Fields 1: " << fields[1] << std::endl;
        std::cerr << "pop1: " << p[fields[0]] << " pop 2: " << p[fields[1]] << std::endl;
#endif
        doubletons[p_int_map[p[fields[0]]]][p_int_map[p[fields[1]]]]++;
    }
    if (result.counts.overall == (numChromosomes - 2)) {
        std::vector<int>::iterator it_hom;
        std::vector<int>::iterator it_het;
        // Check if the doubleton is in a single homozygous individual
        it_hom = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 0);
        if (it_hom != result.counts.individualsWithVariant.end()) {
            int doubleton_hom = std::distance(result.counts.individualsWithVariant.begin(), it_hom);
            fields[0] = doubleton_hom; fields[1] = doubleton_hom;
        } else {
            // Otherwise look at which two individuals have the two heterozygous variants
            it_het = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 1);
            int i = 0;
            while (it_het != result.counts.individualsWithVariant.end()) {
                int doubleton_het = std::distance(result.counts.individualsWithVariant.begin(), it_het);
                fields[i] = doubleton_het;
                ++it_het;
                it_het = std::find(it_het, result.counts.individualsWithVariant.end(), 1);
                ++i;
            }
            assert(i <= 3);
        }
#ifdef TESTING
        std::cerr << "Fields 0: " << fields[0] << " Fields 1: " << fields[1] << std::endl;
        std::cerr << "pop1: " << p[fields[0]] << " pop 2: " << p[fields[1]] << std::endl;
#endif
        doubletons[p_int_map[p[fields[0]]]][p_int_map[p[fields[1]]]]++;
    }
}



// Increment the diff (half)matrix, recording the differences 
void diffs_between_individuals(std::vector<std::vector<double> >& diffs,std::vector<std::vector<double> >& diffs_me, std::vector<std::vector<double> >& diffs_Hets_vs_Homs, FilterResult& result) {
    for (std::vector<std::vector<int> >::size_type i = 0; i < diffs.size(); i++) {
        for (int j = 0; j <= i; j++) {
            const int ind_i = result.counts.individualsWithVariant[i];
            const int ind_j = result.counts.individualsWithVariant[j];
            double diff_measure_Richard; double diff_measure_Me;
            
            
            if (j < i) {
                if (ind_i == 1 && ind_j == 1) { // a pair of individuals are both heterozygous for a variant
                    diff_measure_Richard = DIFF_WEIGHT_BOTH_HETS_RICHARD;
                    diff_measure_Me = DIFF_WEIGHT_BOTH_HETS_ME;
                    diffs_Hets_vs_Homs[i][j]++;  // Cases where both are heterozygous are added in the bottom left corner
                } else if ((ind_i == 2 && ind_j == 0) || (ind_i == 0 && ind_j == 2)) { // one is hom for reference allele the other is hom for alternative allele
                    diff_measure_Richard = DIFF_WEIGHT_HOM_DIFFERENCE_RICHARD;
                    diff_measure_Me = DIFF_WEIGHT_HOM_DIFFERENCE_ME;
                    diffs_Hets_vs_Homs[j][i]++;  // Homozygous differences are added in the top right corner
                }
                // one is homozygous the other is heterozygous for alternative allele
                else if ((ind_i == 2 && ind_j == 1) || (ind_i == 1 && ind_j == 2) || (ind_i == 0 && ind_j == 1) || (ind_i == 1 && ind_j == 0)) {
                    diff_measure_Richard = DIFF_WEIGHT_ONE_HOM_ONE_HET;
                    diff_measure_Me = DIFF_WEIGHT_ONE_HOM_ONE_HET;
                } else {
                    diff_measure_Richard = 0;
                    diff_measure_Me = 0;
                }   
                diffs[i][j] = diffs[i][j] + diff_measure_Richard;
                diffs_me[i][j] = diffs_me[i][j] + diff_measure_Me;
            } else if (j == i) { // Fill in the diagonal, the number of hets
                if (ind_i == 1) {
                    diffs[i][j]++;
                    diffs_me[i][j]++;
                }
            }
        }
    }
}