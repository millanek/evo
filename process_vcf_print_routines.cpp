//
//  process_vcf_print_routines.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_print_routines.h"


// Printing doubletons
void print_doubleton_distribution(const string& fileRoot, const std::vector<std::string>& header, std::vector<std::vector<int> >& doubletons) {
    rearrange_doubletons(doubletons);
    std::ios_base::openmode mode_out = std::ios_base::out;
    string doubletonFileName = fileRoot + ".doubletons.txt";
    std::ofstream* pDoubletonOutFile = new std::ofstream(doubletonFileName.c_str(), mode_out);
    *pDoubletonOutFile << "# Doubleton distribution:" << fileRoot << ".vcf" << std::endl;
    *pDoubletonOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    
    // Print the doubletons matrix
    print_vector(header,*pDoubletonOutFile);
    print_matrix<std::vector<std::vector<int> >&>(doubletons, *pDoubletonOutFile);
    
}

// Printing het counts
void print_het_counts(const string& fileRoot, const std::vector<std::string>& header, const std::vector<int>& hetCounts, const std::vector<int>& sharedHetCounts) {
    assert(hetCounts.size() == sharedHetCounts.size());
    std::ios_base::openmode mode_out = std::ios_base::out;
    string hetFileName = fileRoot + ".hets.txt";
    string sharedHetFileName = fileRoot + ".sharedHets.txt";
    std::ofstream* pHetsOutFile = new std::ofstream(hetFileName.c_str(), mode_out);
    std::ofstream* pSharedHetsOutFile = new std::ofstream(sharedHetFileName.c_str(), mode_out);
    *pHetsOutFile << "# Het counts" << std::endl;
    *pHetsOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    
    *pSharedHetsOutFile << "# Shared het counts (line1) and proportions (line 2)" << std::endl;
    *pSharedHetsOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;

    // Calculate shared het proportions
    std::vector<double> sharedHetProportions;
    for (int i = 0; i < (int)sharedHetCounts.size(); i++) {
        sharedHetProportions.push_back((double)sharedHetCounts[i]/hetCounts[i]);
    }
    
    // print het counts
    print_vector(header,*pHetsOutFile);
    print_vector(hetCounts, *pHetsOutFile);
    
    // print shared het counts and proportions
    print_vector(header,*pSharedHetsOutFile);
    print_vector(sharedHetCounts, *pSharedHetsOutFile);
    print_vector(sharedHetProportions, *pSharedHetsOutFile);
}


void print_privateFixedVarsSummary(const string& fileRoot, const std::vector<std::string>& header, const string& populationsFile,  const std::vector<int>& privateVarCounts) {
    string privateVarFileName = fileRoot + "_" + stripExtension(populationsFile) + ".privateFixedVars.txt";
    std::ofstream* pPrivateVarFile = new std::ofstream(privateVarFileName.c_str());
    *pPrivateVarFile << "# Counts of private fixed variants:" << std::endl;
    *pPrivateVarFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    *pPrivateVarFile << "# Groups defined in:" << populationsFile << std::endl;
    
    // print het counts
    print_vector(header, *pPrivateVarFile);
    print_vector(privateVarCounts, *pPrivateVarFile);
}



// Printing pairwise difference statistics
void print_pairwise_diff_stats(const string& fileRoot, const std::vector<std::string>& header, const int totalVariantNumber, const std::vector<std::vector<double> >& diffMatrix, const std::vector<std::vector<double> >& diffMatrixMe, const std::vector<std::vector<double> >& diffMatrixHetsVsHomDiff) {
    std::ios_base::openmode mode_out = std::ios_base::out;
    string diffFileName = fileRoot + ".diff_matrix.txt";
    string diffMeFileName = fileRoot + ".diff_me_matrix.txt";
    string hetHomFileName = fileRoot + ".hets_over_homs_matrix.txt";
    std::ofstream* pDiffOutFile = new std::ofstream(diffFileName.c_str(), mode_out);
    std::ofstream* pDiffMeOutFile = new std::ofstream(diffMeFileName.c_str(), mode_out);
    std::ofstream* pHetHomOutFile = new std::ofstream(hetHomFileName.c_str(), mode_out);
    *pDiffOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    *pDiffOutFile << "# Total number of segragating variant sites in this sample:" << totalVariantNumber << std::endl;
    *pDiffOutFile << "# Richard's scoring scheme" << std::endl;
    *pDiffMeOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    *pDiffMeOutFile << "# Total number of segragating variant sites in this sample: " << totalVariantNumber << std::endl;
    *pDiffMeOutFile << "# Homozygous difference = 2, one homozygous, another heterozygous = 1:" << totalVariantNumber << std::endl;
    *pHetHomOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    *pHetHomOutFile << "# number of sites both individuals hets/number of sites individuals have a homozygous difference; i.e. num(1/0::1/0)/num(1/1::0/0)" << std::endl;
    *pHetHomOutFile << "# For a free mixing population, we expect this number ~2; for fully separated species ~0" << std::endl;
    
    // print headers
    print_vector(header,*pDiffOutFile);
    print_vector(header,*pDiffMeOutFile);
    print_vector(header,*pHetHomOutFile);
    
    // print statistics
    print_matrix<const std::vector<std::vector<double> >&>(diffMatrix, *pDiffOutFile);
    print_matrix<const std::vector<std::vector<double> >&>(diffMatrixMe, *pDiffMeOutFile);
    print_matrix<const std::vector<std::vector<double> >&>(diffMatrixHetsVsHomDiff, *pHetHomOutFile);
    
}

// Printing haplotype pairwise difference statistics
void print_H1_pairwise_diff_stats(const string& fileRoot, std::vector<std::string>& header, const int totalVariantNumber, const std::vector<std::vector<double> >& diffMatrixH1) {
    std::ios_base::openmode mode_out = std::ios_base::out;
    string diffFileNameH1 = fileRoot + ".diff_matrix_H1.txt";
    std::ofstream* pDiffH1OutFile = new std::ofstream(diffFileNameH1.c_str(), mode_out);
    *pDiffH1OutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    *pDiffH1OutFile << "# Total number of segragating variant sites in this sample:" << totalVariantNumber << std::endl;
    *pDiffH1OutFile << "# Differences between H1 haplotypes:" << std::endl;

    for (std::vector<std::string>::size_type i = 0; i < header.size(); i++) {
        header[i] = header[i] + "_H1";
    }
    
    // print
    print_vector(header,*pDiffH1OutFile);
    print_matrix<const std::vector<std::vector<double> >&>(diffMatrixH1, *pDiffH1OutFile);
    
}

// Printing haplotype pairwise difference statistics
void print_AllH_pairwise_diff_stats(const string& fileRoot, const std::vector<std::string>& samples, const int totalVariantNumber, const std::vector<std::vector<double> >& diffMatrixAllH) {
    std::ios_base::openmode mode_out = std::ios_base::out;
    string diffFileNameAllH = fileRoot + ".diff_matrix_AllH.txt";
    std::ofstream* pDiffAllHOutFile = new std::ofstream(diffFileNameAllH.c_str(), mode_out);
    *pDiffAllHOutFile << "# Input file:" << fileRoot << ".vcf" << std::endl;
    *pDiffAllHOutFile << "# Total number of segragating variant sites in this sample:" << totalVariantNumber << std::endl;
    *pDiffAllHOutFile << "# Differences between all haplotypes:" << std::endl;
    
    std::vector<std::string> header;
    for (std::vector<std::string>::size_type i = 0; i < samples.size(); i++) {
        header.push_back(samples[i] + "_H1");
        header.push_back(samples[i] + "_H2");
    }
    
    // print statistics
    print_vector(header,*pDiffAllHOutFile);
    print_matrix<const std::vector<std::vector<double> >&>(diffMatrixAllH, *pDiffAllHOutFile);
    
}



