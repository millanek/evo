//
//  extract_target_sequence_from_UTR.h
//  vcf_process
//
//  Created by Milan Malinsky on 20/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef vcf_process_extract_target_sequence_from_UTR_h
#define vcf_process_extract_target_sequence_from_UTR_h

void parseUTROptions(int argc, char** argv);


// ------------------------- FUNCTIONS -------------------------------------
int countSegregatingSites(const std::vector<std::string>& utrVector);
void processAnnotationFile(const std::string& line, std::map<std::string,std::vector<std::string> >& annotationMap);
int writeMzebraCoordinatesIntoBedFiles(const std::string& fileRoot, std::map<std::string,std::vector<std::string> >& annotationMap, const std::vector<int>& targetLoci, std::vector<int>::size_type i, std::ofstream*& targetBed, std::ofstream*& nonTargetBed, int lettersBeforeCutStart, int targetLength, int nextBaseBed);
void writeLastMzebraCoordinateToBed(const std::string& fileRoot, std::map<std::string,std::vector<std::string> >& annotationMap, const std::vector<int>& targetLoci, std::ofstream*& nonTargetBed, int lettersBeforeCutStart, int nextBaseBed, std::string::size_type l);
// ------------------------- UTILS -------------------------------------
// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename);
// Split strings
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
#endif
