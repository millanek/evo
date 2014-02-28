//
//  filter_annotation.h
//  vcf_process
//
//  Created by Milan Malinsky on 22/07/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef vcf_process_filter_annotation_h
#define vcf_process_filter_annotation_h

// ------------------------- UTILS -------------------------------------
// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename);
// Split strings
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);


#endif
