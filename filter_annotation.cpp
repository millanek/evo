//
//  filter_annotation.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 22/07/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <map>
#include "filter_annotation.h"
using std::string;

int main(int argc, char **argv) {
    std::map<string,string> brawandMap;
    string brawandFileName = "/Users/milanmalinsky/Work/annotation/ci_orthologous_clusters.txt";
    std::ifstream* brawandFile = new std::ifstream(brawandFileName.c_str()); 
    string line;
    while (getline(*brawandFile, line)) {
        if (line[0] == 'm' && line[1] == 'z') {
            std::vector<string> nameIdVec = split(line, '\t');
            brawandMap[nameIdVec[0]] = nameIdVec[1];
        }
    }
    
    string annotationFileName = "/Users/milanmalinsky/Work/annotation/annotation_files/Metriaclima_zebra.BROADMZ1.gff3";
    std::ifstream* annotationFile = new std::ifstream(annotationFileName.c_str()); 
    while (getline(*annotationFile, line)) {
        std::vector<string> annotLineVec = split(line, '\t');
        std::vector<string> lineDescVec = split(annotLineVec[8], ';');
        std::vector<string> geneNameVec = split(lineDescVec[lineDescVec.size()-1], '.');
        string thisGene = "mz.gene." + geneNameVec[2] + "." + geneNameVec[3];
        if (brawandMap.count(thisGene) != 0) {
            std::cout << line << std::endl;
        }
    }
}



// ------------------------- UTILS -------------------------------------

// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}


void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
