//
//  process_vcf_stats_functions.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_stats_functions.h"

static const double DIFF_WEIGHT_BOTH_HETS_ME = 0.5;
static const double DIFF_WEIGHT_BOTH_HETS_RICHARD = (2.0/3.0);
static const double DIFF_WEIGHT_HOM_DIFFERENCE_ME = 1;
static const double DIFF_WEIGHT_HOM_DIFFERENCE_RICHARD = (2.0/3.0);
static const double DIFF_WEIGHT_ONE_HOM_ONE_HET = 0.5;

// Increment a counter for each individual who is het at this site
void het_analysis(std::vector<int>& hetCounts, std::vector<int>& sharedHetCounts, FilterResult& result) {
    for (std::vector<int>::size_type i = 0; i < hetCounts.size(); i++) {
        if (result.counts.individualsWithVariant[i] == 1) {
            hetCounts[i]++;
            if (result.counts.overall > 1) {
                sharedHetCounts[i]++;
            }
        }
    }
}

std::vector<string> initialisePopulationMap(const std::vector<std::string> populations, std::map<std::string,int>& fp_map) {
    std::vector<std::string> pop_unique = populations;
    std::sort(pop_unique.begin(), pop_unique.end());
    std::vector<std::string>::iterator it = std::unique(pop_unique.begin(), pop_unique.end());
    pop_unique.resize(std::distance(pop_unique.begin(), it));
    std::vector<std::string>::size_type n_pop = pop_unique.size();
    // Initialize fields to populations map
    for (int i = 0; i < n_pop; i++) {
        fp_map[pop_unique[i]] = i;
    }
    std::cerr << "Initialised the population map; there are " << n_pop << " populations" << std::endl;
    return pop_unique;
}

std::map<std::string,std::vector<int>> getPopulationsToIndividualsMap(std::vector<std::string> populations, std::map<std::string,int>& fp_map) {
    std::map<std::string,std::vector<int>> popsToIndivMap;
    for (std::map<std::string,int>::iterator it = fp_map.begin(); it != fp_map.end(); it++) {
        std::vector<int> thisPopI;
        std::vector<std::string>::iterator it_i =  std::find(populations.begin(), populations.end(), it->first);
        while (it_i != populations.end()) {
            thisPopI.push_back((int)std::distance(populations.begin(), it_i));
            ++it_i;
            it_i = std::find(it_i, populations.end(), it->first);
        }
        popsToIndivMap[it->first] = thisPopI;
    }
    return popsToIndivMap;
}




// Initializing the doubleton data structure
void initializeDoubletons(std::vector<std::vector<int> >& d, const std::vector<std::string> populationsUnique) {
    std::vector<std::string>::size_type n_pop = populationsUnique.size();
    
    // Initialize doubletons
    for (int i = 0; i < n_pop; i++) { 
        std::vector<int> v(n_pop,0);
        d.push_back(v);
    }
#ifdef TESTING
    std::cerr << "Doubletons initialised with size: " << n_pop << std::endl;
#endif
}


void privateVars_analysis(std::vector<int>& privateVarCounts, const FilterResult& result, const std::vector<std::vector<size_t> >& populationsIndices, const std::vector<std::vector<size_t> >& populationsIndicesComplements) {
    for (int i = 0; i < (int)populationsIndices.size(); i++) {
        bool allAlts = true; bool allRefs = true;
        for (int j = 0; j < populationsIndices[i].size(); j++) {
            if (result.counts.individualsWithVariant[populationsIndices[i][j]] != 2)
                allAlts = false;
            if (result.counts.individualsWithVariant[populationsIndices[i][j]] != 0)
                allRefs = false;
            if (!allAlts && !allRefs)
                break;
        }
        if (!allAlts && !allRefs)
            break;
        
        // If all individuals in the population have the same allele (either all alt or all ref)
        // Then we look at the rest of the individuals in the dataset (complement)
        bool privateVar = true;
        if (allAlts) {
            for (int j = 0; j < populationsIndicesComplements[i].size(); j++) {
                if (result.counts.individualsWithVariant[populationsIndicesComplements[i][j]] != 0)
                    privateVar = false;
            }
            //print_vector_stream(populationsIndices[i], std::cerr);
            //print_vector_stream(result.counts.individualsWithVariant, std::cerr);
        }
        if (allRefs) {
            for (int j = 0; j < (int)populationsIndicesComplements[i].size(); j++) {
                if (result.counts.individualsWithVariant[populationsIndicesComplements[i][j]] != 2)
                privateVar = false;
            }
        }
        if (privateVar) {
            privateVarCounts[i]++;
        }
    }
}

// Filling in the doubleton data structure
void doubleton_analysis(std::vector<std::vector<int> >& doubletons, FilterResult& result, int numChromosomes, const std::vector<std::string> p, std::map<std::string,int>& p_int_map) {
#ifdef TESTING
    std::cerr << "Overall count: " << result.counts.overall << std::endl;
#endif
    int fields[2];
    
    if (result.counts.overall == 2) {
        std::vector<int>::iterator it_hom;
        std::vector<int>::iterator it_het;
        // Check if the doubleton is in a single homozygous individual
        it_hom = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 2);
        if (it_hom != result.counts.individualsWithVariant.end()) {
            int doubleton_hom = (int)std::distance(result.counts.individualsWithVariant.begin(), it_hom);
            fields[0] = doubleton_hom; fields[1] = doubleton_hom;
        } else {
            // Otherwise look at which two individuals have the two heterozygous variants
            it_het = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 1);
            int i = 0;
            while (it_het != result.counts.individualsWithVariant.end()) {
                int doubleton_het = (int)std::distance(result.counts.individualsWithVariant.begin(), it_het);
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
            int doubleton_hom = (int)std::distance(result.counts.individualsWithVariant.begin(), it_hom);
            fields[0] = doubleton_hom; fields[1] = doubleton_hom;
        } else {
            // Otherwise look at which two individuals have the two heterozygous variants
            it_het = std::find(result.counts.individualsWithVariant.begin(), result.counts.individualsWithVariant.end(), 1);
            int i = 0;
            while (it_het != result.counts.individualsWithVariant.end()) {
                int doubleton_het = (int)std::distance(result.counts.individualsWithVariant.begin(), it_het);
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
void diffs_between_individuals(std::vector<std::vector<double> >& diffs,std::vector<std::vector<double> >& diffs_me, std::vector<std::vector<double> >& diffs_me_bootstrap, std::vector<std::vector<double> >& diffs_Hets_vs_Homs, std::vector<std::vector<int> >& pairwise_missingness, std::vector<std::vector<int> >& pairwise_missingness_bootstrap, FilterResult& result) {
    for (std::vector<std::vector<int> >::size_type i = 0; i < diffs.size(); i++) {
        const int ind_i = result.counts.individualsWithVariant[i];
        const bool ind_i_missing = result.counts.missingGenotypesPerIndividual[i];
        for (int j = 0; j <= i; j++) {
            double diff_measure_Richard; double diff_measure_Me;
            if (ind_i_missing) {
                diff_measure_Me = 0;
                diff_measure_Richard = 0;
                pairwise_missingness[i][j]++;
                pairwise_missingness_bootstrap[i][j]++;
            } else {
                const bool ind_j_missing = result.counts.missingGenotypesPerIndividual[j];
                if (ind_j_missing) {
                    diff_measure_Me = 0;
                    diff_measure_Richard = 0;
                    pairwise_missingness[i][j]++;
                    pairwise_missingness_bootstrap[i][j]++;
                } else {
                    const int ind_j = result.counts.individualsWithVariant[j];
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
                        diffs_me_bootstrap[i][j] = diffs_me_bootstrap[i][j] + diff_measure_Me;
                    } else if (j == i) { // Fill in the diagonal, the number of hets
                        if (ind_i == 1) {
                            diffs[i][j]++;
                            diffs_me[i][j]++;
                            diffs_me_bootstrap[i][j]++;
                        }
                    }
                }
            }
        }
    }
}

// Increment the diff (half)matrix, recording the differences
void diffs_between_individuals_with_multialleleics(std::vector<std::vector<double> >& diffs_me, std::vector<std::vector<int> >& pairwise_missingness, std::vector<std::vector<double> >& diffs_me_bootstrap, std::vector<std::vector<int> >& pairwise_missingness_bootstrap, FilterResult& result) {
    for (std::vector<std::vector<int> >::size_type i = 0; i < diffs_me.size(); i++) {
        const int ind_i = result.counts.haplotypesWithVariant[2*i];
        const int ind_i2 = result.counts.haplotypesWithVariant[2*i+1];
        const bool ind_i_missing = result.counts.missingGenotypesPerIndividual[i];
        for (int j = 0; j <= i; j++) {
            double diff_measure_Me;
            if (ind_i_missing) {
                diff_measure_Me = 0;
                pairwise_missingness[i][j]++;
                pairwise_missingness_bootstrap[i][j]++;
            } else {
                const bool ind_j_missing = result.counts.missingGenotypesPerIndividual[j];
                if (ind_j_missing) {
                    diff_measure_Me = 0;
                    pairwise_missingness[i][j]++;
                    pairwise_missingness_bootstrap[i][j]++;
                } else {
                    const int ind_j = result.counts.haplotypesWithVariant[2*j];
                    const int ind_j2 = result.counts.haplotypesWithVariant[2*j+1];
                    if (j < i) {
                        double totalD = 0; const double numComparisons = 4;
                        if (ind_i != ind_j) totalD++;
                        if (ind_i != ind_j2) totalD++;
                        if (ind_i2 != ind_j) totalD++;
                        if (ind_i2 != ind_j2) totalD++;
                        diff_measure_Me = totalD/numComparisons;
                        if (diff_measure_Me != 0 && diff_measure_Me != 0.5 && diff_measure_Me != 0.75 && diff_measure_Me != 1) {
                            std::cerr << "diff_measure_Me: " << diff_measure_Me << std::endl;
                            std::cerr << "ind_i: " << ind_i << std::endl;
                            std::cerr << "ind_i2: " << ind_i2 << std::endl;
                            std::cerr << "ind_j: " << ind_j << std::endl;
                            std::cerr << "ind_j2: " << ind_j2 << std::endl;
                        }
                        diffs_me[i][j] = diffs_me[i][j] + diff_measure_Me;
                        diffs_me_bootstrap[i][j] = diffs_me_bootstrap[i][j] + diff_measure_Me;
                    } else if (j == i) { // Fill in the diagonal, the number of hets
                        if (ind_i != ind_i2) {
                            diffs_me[i][j]++;
                            diffs_me_bootstrap[i][j]++;
                        }
                    }
                }
            }
        }
    }
}



void diffs_between_H1(std::vector<std::vector<double> >& H1diffs, FilterResult& result) {
    for (std::vector<std::vector<int> >::size_type i = 0; i < H1diffs.size(); i++) {
        for (int j = 0; j < i; j++) {
            const int ind_i = result.counts.haplotypesWithVariant[2*i];
            const int ind_j = result.counts.haplotypesWithVariant[2*j];
            
            if ((ind_i == 1 && ind_j == 0) || (ind_i == 0 && ind_j == 1)) {
                H1diffs[i][j]++;
            } else if ((ind_i == 0 && ind_j == 0) || (ind_i == 1 && ind_j == 1)) { // one is hom for reference allele the other is hom for alternative allele
                
            } else {
                std::cerr << "result.counts.haplotypesWithVariant[" << 2*i << "]=" << result.counts.haplotypesWithVariant[2*i] << std::endl;
                std::cerr << "result.counts.haplotypesWithVariant[" << 2*j << "]=" << result.counts.haplotypesWithVariant[2*j] << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}

void diffs_between_AllH(std::vector<std::vector<double> >& AllHdiffs, FilterResult& result) {
    for (std::vector<std::vector<int> >::size_type i = 0; i < AllHdiffs.size(); i++) {
        for (int j = 0; j < i; j++) {
            const int ind_i = result.counts.haplotypesWithVariant[i];
            const int ind_j = result.counts.haplotypesWithVariant[j];
            
            if ((ind_i == 1 && ind_j == 0) || (ind_i == 0 && ind_j == 1)) {
                AllHdiffs[i][j]++;
            } else if ((ind_i == 0 && ind_j == 0) || (ind_i == 1 && ind_j == 1)) { // one is hom for reference allele the other is hom for alternative allele
                
            } else {
                std::cerr << "result.counts.haplotypesWithVariant[" << i << "]=" << result.counts.haplotypesWithVariant[i] << std::endl;
                std::cerr << "result.counts.haplotypesWithVariant[" << j << "]=" << result.counts.haplotypesWithVariant[j] << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}

