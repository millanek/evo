//
//  process_vcf_sequenom.h
//  vcf_process
//
//  Created by Milan Malinsky on 29/08/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef vcf_process_process_vcf_sequenom_h
#define vcf_process_process_vcf_sequenom_h


void writeSequenomOutput(const std::string& regionSeq, const std::string& refSeq, const std::string& loc, std::ofstream*& sequenomFile);

int sequenomMain(int argc, char** argv);
void parseSequenomOptions(int argc, char** argv);

#endif
