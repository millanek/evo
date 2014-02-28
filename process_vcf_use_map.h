//
//  process_vcf_use_map.h
//  vcf_process
//
//  Created by Milan Malinsky on 08/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_use_map__
#define __vcf_process__process_vcf_use_map__

#include <iostream>


std::map<std::string, std::vector<std::vector<std::string> > > loadLinkageGroupMap(std::ifstream*& LGfile);

void parseMapOptions(int argc, char** argv);

int mapMain(int argc, char** argv);
void processGenome(const std::map<std::string, std::vector<std::vector<std::string> > >& LGmap);
void processVCF(const std::map<std::string, std::vector<std::vector<std::string> > >& LGmap);
void printProcessedVCFline(std::vector<string>& vcfFields, const bool bInLG, const char thisScaffoldOrientation, const std::string& thisLG, const int& LGsizeUpToHere, const int& thisScaffoldSize, const std::string& line);

#endif /* defined(__vcf_process__process_vcf_use_map__) */
