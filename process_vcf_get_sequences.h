//
//  process_vcf_get_sequences.h
//  vcf_process
//
//  Created by Milan Malinsky on 17/07/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef vcf_process_process_vcf_get_sequences_h
#define vcf_process_process_vcf_get_sequences_h

#include "process_vcf_seq_utils.h"

int getSeqMain(int argc, char** argv);
void parseGetSeqOptions(int argc, char** argv);
void printInAllOutputs(std::ofstream*& outFiles, size_t numSamples, std::string toPrint);
void print80bpPerLine(std::ofstream*& outFiles, std::vector<std::string>::size_type i, std::string toPrint);
void print_split(const std::string& currentScaffoldNum, const std::vector<string::size_type>& splits, const std::vector<std::string>& sampleNames, const size_t numSamples, std::vector<std::string>& scaffoldStrings, const unsigned int totalProcessedVariants);
void print_split_incl_outgroup(const std::string& currentScaffoldNum, const std::vector<string::size_type>& splits, const std::vector<std::string>& sampleNames, const size_t numSamples, std::vector<std::string>& scaffoldStrings, const unsigned int totalProcessedVariants, std::map<std::string, std::string>& outgroupSeqs, const std::string& outgroupName);

#endif
