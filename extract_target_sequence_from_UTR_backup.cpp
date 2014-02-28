//
//  extract_target_sequence_from_UTR.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 20/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <getopt.h>
#include <cstdlib>
#include <dirent.h>
#include "extract_target_sequence_from_UTR.h"
using std::string;
#define PROGRAM_BIN "extract_target_sequence_from_UTR"
#define PACKAGE_BUGREPORT "mm812@cam.ac.uk"

//#define TESTING 1


#define AUTHOR "Milan Malinsky"
#define PACKAGE_VERSION "0.1"


static const char *USAGE_MESSAGE =
"Program: " PROGRAM_BIN "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"Commands:\n"
"           filter              filter a vcf file\n"
"           stats               get various statistics from a vcf file\n"
"           sharing             find out how Massoko variants segregate in Lake Malawi\n"
"           test                used for testing/developing new programs\n"
"           massoko             filtering and getting fixed (and nearly fixed) variants from the Lake Massoko VCF file"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


namespace opt
{
    static string utrFile;
    static string targetLociFile;
}

int main(int argc, char **argv) {
    
    parseUTROptions(argc, argv);
    
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir ("/Users/milanmalinsky/Work/annotation/multiple_align_UTRs")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            string thisFileName(ent->d_name);
            std::vector<string> fname_parts = split(thisFileName, '.');
            if (fname_parts[fname_parts.size()-1] == "dashrem") {
                printf ("%s\n", ent->d_name);
            }
        }
        closedir (dir);
    } else {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }
    
    std::ios_base::openmode mode = std::ios_base::in;
    string utrFileName = opt::utrFile;
    string fileRoot = stripExtension(utrFileName);
    string targetFileName = opt::targetLociFile;
    std::vector<string> UTRs;
    std::vector<int> targetLoci;
    std::map<int,int> mzTranslateCoordinates;
    
    // Open connection to read from the utr and target loci files
    std::ifstream* utrInFile = new std::ifstream(utrFileName.c_str(), mode);
    std::ifstream* targetInFile = new std::ifstream(targetFileName.c_str(), mode);
    
    
    string line;
    while (getline(*utrInFile, line)) {
        std::vector<string> species_UTR = split(line, '\t');
        UTRs.push_back(species_UTR[1]);
        // std::cout << line << std::endl;
    }
    
    while (getline(*targetInFile, line)) {
        targetLoci.push_back(atoi(line.c_str())-1);
        // std::cout << line << std::endl;
    }
    // Get unique target loci
    std::sort(targetLoci.begin(), targetLoci.end());
    std::vector<int>::iterator it = std::unique (targetLoci.begin(), targetLoci.end());
    targetLoci.resize(std::distance(targetLoci.begin(), it));
    
    string mzeb_UTR = UTRs[0];
    
    // Prepare a map to have translating from UTR coordinates to aligned (with dashes==indels)
    int dashes_so_far = 0;
    int letters_so_far = 0;
    for (string::size_type i = 0; i != mzeb_UTR.length(); i++) {
        if (mzeb_UTR[i] == '-')
            dashes_so_far++;
        mzTranslateCoordinates[letters_so_far] = letters_so_far + dashes_so_far;
        if (mzeb_UTR[i] == 'A' || mzeb_UTR[i] == 'C' || mzeb_UTR[i] == 'G' || mzeb_UTR[i] == 'T')
            letters_so_far++;
    }
    
   // for (std::map<int, int>::iterator i = mzTranslateCoordinates.begin(); i != mzTranslateCoordinates.end(); i++) {
    //    std::cout << "Unaligned pos: " << i->first << " Aligned pos: " << i->second << std::endl;
    //}
    
    for (std::vector<string>::size_type j = 0; j != UTRs.size(); j++) {
        for (std::vector<int>::size_type i = 0; i != targetLoci.size(); i++) {
            if (mzTranslateCoordinates.find(targetLoci[i]) == mzTranslateCoordinates.end()) {
                // std::cout << "Target " << i+1 << " without dashes: " << targetLoci[i] << " is beyond the alignment region" << std::endl; 
            } else {
                // std::cout << "Target " << i+1 << " without dashes: " << targetLoci[i] << " Locus start: " << mzTranslateCoordinates[targetLoci[i]] << std::endl;  
                std::cout << UTRs[j].substr(mzTranslateCoordinates[targetLoci[i]],7);
            }
        }
        if (targetLoci.size() > 0) 
            std::cout << std::endl;
    }
    return 0;
    
}


void parseUTROptions(int argc, char** argv) {
    bool die = false;
   /* for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'd': arg >> opt::max_overall_depth; break;
            case 'c': arg >> opt::min_copies; break;
            case 's': arg >> opt::min_depth_in_any_individual; break;
            case '?': die = true; break;   
            case OPT_HELP:
                std::cout << FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    } */
      
    if (argc != 3) 
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::utrFile = argv[1];
    opt::targetLociFile = argv[2];
    
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
