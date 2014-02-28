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
// #define DEBUG1 1
// #define DEBUG2 1

static const char *USAGE_MESSAGE =
"Program: " PROGRAM_BIN "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


namespace opt
{
    static string utrFile;
    static string targetLociFile;
}

int main(int argc, char **argv) {
    
    parseUTROptions(argc, argv);
    std::ios_base::openmode mode = std::ios_base::in;
    
    // Read in the M. zebra annotation file
    std::map<string,std::vector<string> > annotationMap;
    string annotationFileName = "/Users/milanmalinsky/Work/annotation/present_in_all_cichlids_3UTR_cichlids_updated_all_sorted_fullCollapsed.fasta";
    std::ifstream* annotationFile = new std::ifstream(annotationFileName.c_str(), mode); 
    string line;
    while (getline(*annotationFile, line)) {
        if (line[0] == '>') {
            processAnnotationFile(line,annotationMap);
        }
    }
    
    
    // Open connections to write bed files with target and non-target loci
    string targetBedName = "target_loci.bed";
    string targetAlignmentName = "target_alignment.aln";
    string nonTargetBedName = "nonTarget_loci.bed";
    string nonTargetAlignmentName = "nonTarget_alignment.aln";
    std::ofstream* targetBed = new std::ofstream(targetBedName.c_str());
    std::ofstream* nonTargetBed = new std::ofstream(nonTargetBedName.c_str());
    std::ofstream* targetAlignment = new std::ofstream(targetAlignmentName.c_str());
    std::ofstream* nonTargetAlignment = new std::ofstream(nonTargetAlignmentName.c_str());
    
    // Initialise vectors to hold concatenated target and non-target sequences
    std::vector<string> targetsConcat;
    std::vector<string> nonTargetConcat;
    string mzTargetsConcat = ""; targetsConcat.push_back(mzTargetsConcat); nonTargetConcat.push_back(mzTargetsConcat);
    string pnTargetsConcat = ""; targetsConcat.push_back(pnTargetsConcat); nonTargetConcat.push_back(pnTargetsConcat);
    string abTargetsConcat = ""; targetsConcat.push_back(abTargetsConcat); nonTargetConcat.push_back(abTargetsConcat);
    string nbTargetsConcat = ""; targetsConcat.push_back(nbTargetsConcat); nonTargetConcat.push_back(nbTargetsConcat);
    string onTargetsConcat = ""; targetsConcat.push_back(onTargetsConcat); nonTargetConcat.push_back(onTargetsConcat);
    string starsTargetsConcat = ""; targetsConcat.push_back(starsTargetsConcat); nonTargetConcat.push_back(starsTargetsConcat);
    
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir ("/Users/milanmalinsky/Work/annotation/multiple_align_UTRs")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            string thisFileName(ent->d_name);
            std::vector<string> fname_parts = split(thisFileName, '.');
            if (fname_parts[fname_parts.size()-1] == "dashrem") {
                string utrFileName = thisFileName;
                string fileRoot = stripExtension(utrFileName);
                string targetFileName = fileRoot + ".targetloci";
                string fileStartCut = fileRoot + ".dashlstartmax";
                std::vector<string> UTRs;
                std::vector<int> targetLoci;
                std::map<int,int> mzTranslateCoordinates;
                int cutStart = 0;
                
                
                // Open connection to read from the utr and target loci files
                std::ifstream* utrInFile = new std::ifstream(utrFileName.c_str(), mode);
                std::ifstream* targetInFile = new std::ifstream(targetFileName.c_str(), mode);
                std::ifstream* cutInFile = new std::ifstream(fileStartCut.c_str(), mode);
                
                // Read in data
                string line;
                while (getline(*utrInFile, line)) {
                    std::vector<string> species_UTR = split(line, '\t');
                    UTRs.push_back(species_UTR[1]);
                   // std::cout << line << std::endl;
                }
                
                getline(*cutInFile, line);
                cutStart = atoi(line.c_str());
                // std::cout << cutStart << std::endl;
                // Ignore alignments where, after removing overhangs from the beginning and the end, nothing is left
                if (cutStart > UTRs[0].length()) { // In these cases there is no sequence left
                    std::cout << "Error in: " << thisFileName << std::endl;
                    std::cout << "CutStart: " << cutStart << " UTR length: " << UTRs[0].length() << std::endl;
                    continue;
                }
                
                
                while (getline(*targetInFile, line)) {
                    targetLoci.push_back(atoi(line.c_str())-1);
                    //std::cout << line << std::endl;
                }
                // Get unique target loci
                std::sort(targetLoci.begin(), targetLoci.end());
                std::vector<int>::iterator it = std::unique (targetLoci.begin(), targetLoci.end());
                targetLoci.resize(std::distance(targetLoci.begin(), it));
                
                
                // Sometimes, two targets are closer than 7bp from each other (i.e. they overlap)
                std::map<int,int> locusLength;
                if (targetLoci.size() > 0) {
                    if (targetLoci.size() == 1) {
                        locusLength[targetLoci[0]] = 7;
                    } else {
                        for (std::vector<int>::size_type i = 1; i != targetLoci.size(); i++) {
                            int diff = targetLoci[i] - targetLoci[i-1];
                            if (diff < 7) 
                                locusLength[targetLoci[i-1]] = diff;
                            else
                                locusLength[targetLoci[i-1]] = 7;
                            if (i == (targetLoci.size()-1)) 
                                locusLength[targetLoci[i]] = 7;
                        }
                    }
                }
                
                // Prepare a map to have translating from UTR coordinates to aligned (with dashes==indels)
                string mzeb_UTR = UTRs[0];
                int dashes_so_far = 0;
                int letters_so_far = 0;
                int lettersBeforeCutStart;
                for (string::size_type i = 0; i != mzeb_UTR.length(); i++) {
                    if (i == cutStart) {
                        lettersBeforeCutStart = letters_so_far;
                    }
                    if (mzeb_UTR[i] == '-')
                        dashes_so_far++;
                    mzTranslateCoordinates[letters_so_far] = letters_so_far + dashes_so_far;
                    if (mzeb_UTR[i] == 'A' || mzeb_UTR[i] == 'C' || mzeb_UTR[i] == 'G' || mzeb_UTR[i] == 'T' || mzeb_UTR[i] == 'N')
                        letters_so_far++;
                }
                
                 
                // Ignore cases 
                
               // for (std::map<int, int>::iterator i = mzTranslateCoordinates.begin(); i != mzTranslateCoordinates.end(); i++) {
                //    std::cout << "Unaligned pos: " << i->first << " Aligned pos: " << i->second << std::endl;
                //}
                
                // Now get the target site sequences
                bool bTargetInAlignment = false;
                for (std::vector<string>::size_type j = 0; j != UTRs.size(); j++) {
                    string::size_type nextBase = cutStart;
                    int nextBaseBed = 0;
                    for (std::vector<int>::size_type i = 0; i != targetLoci.size(); i++) { // For each target in this 3'UTR
                        if (mzTranslateCoordinates.count(targetLoci[i]) == 0) {
                             // std::cout << "Target " << i+1 << " without dashes: " << targetLoci[i] << " is beyond the alignment region" << std::endl;
                        } else {
                            #ifdef DEBUG2
                            std::cout << "Target " << i+1 << " without dashes: " << targetLoci[i] << " Locus start: " << mzTranslateCoordinates[targetLoci[i]] << std::endl;  
                            #endif
                            if (mzTranslateCoordinates[targetLoci[i]] > cutStart) {
                                int targetLength = locusLength[targetLoci[i]];
                                // If this target is entirely within the alignment region
                                if (UTRs[j].length() > mzTranslateCoordinates[targetLoci[i]]+targetLength) {
                                    #ifdef DEBUG1
                                    std::cout << UTRs[j].substr(mzTranslateCoordinates[targetLoci[i]],targetLength);
                                    #endif
                                    targetsConcat[j] = targetsConcat[j] + UTRs[j].substr(mzTranslateCoordinates[targetLoci[i]],targetLength);
                                    nonTargetConcat[j] = nonTargetConcat[j] + UTRs[j].substr(nextBase,mzTranslateCoordinates[targetLoci[i]] - nextBase);
                                    nextBase = mzTranslateCoordinates[targetLoci[i]] + targetLength;
                                    
                                    if (j == 0)  // If this is the M. zebra 3'UTR
                                        nextBaseBed = writeMzebraCoordinatesIntoBedFiles(fileRoot, annotationMap, targetLoci, i, targetBed, nonTargetBed, lettersBeforeCutStart, targetLength, nextBaseBed);
                                    
                                } else {
                                   // std::cout << UTRs[j].substr(mzTranslateCoordinates[targetLoci[i]]);
                                   // targetsConcat[j] = targetsConcat[j] + UTRs[j].substr(mzTranslateCoordinates[targetLoci[i]]);
                                } 
                            } else {
                                // std::cout << "Target " << i+1 << " without dashes: " << targetLoci[i] << " is beyond the alignment region" << std::endl;
                            }
                            bTargetInAlignment = true;
                        }
                    }
                    
                    // Get the reminder of the sequence after the last target in this 3'UTR
                    nonTargetConcat[j] = nonTargetConcat[j] + UTRs[j].substr(nextBase);
                    if (j == 0) {
                        std::string::size_type l = 0;
                        for (string::size_type i = 0; i != UTRs[j].substr(nextBase).length(); i++) {
                            if (UTRs[j].substr(nextBase)[i] != '-') {
                                l++;
                            }
                        }
                        writeLastMzebraCoordinateToBed(fileRoot, annotationMap, targetLoci, nonTargetBed, lettersBeforeCutStart, nextBaseBed, l);
                    }
                    
                    
                    #ifdef DEBUG1
                    if (bTargetInAlignment) 
                        std::cout << std::endl;
                    #endif
                }
                #ifdef DEBUG1
                if (bTargetInAlignment)
                    std::cout << thisFileName << std::endl;
                #endif
                
                // cleanup
                utrInFile->close();
                targetInFile->close();
                cutInFile->close();
                UTRs.clear();
                targetLoci.clear();
                mzTranslateCoordinates.clear();
            }
        }
        closedir (dir);
        
        // Output alignments to the file
        for (std::vector<string>::size_type i = 0; i != targetsConcat.size(); i++) {
            std::cout << targetsConcat[i] << std::endl;
            *targetAlignment << targetsConcat[i] << std::endl;
        }
        for (std::vector<string>::size_type i = 0; i != nonTargetConcat.size(); i++) {
            *nonTargetAlignment << nonTargetConcat[i] << std::endl;
        }
        
        #ifdef DEBUG1
        for (std::map<string, std::vector<string> >::iterator i = annotationMap.begin(); i != annotationMap.end(); i++) {
            std::cout << "Gene: " << i->first << " Locus: " << i->second[0] << "\t" << i->second[1] << "\t" << i->second[2] << std::endl;
        }
        #endif
        int num_same = 0;
        int num_segregating = countSegregatingSites(targetsConcat);
        num_same = targetsConcat[0].length() - num_segregating;
        std::cout << "Alignment in targets length: " << targetsConcat[0].length() << std::endl;
        std::cout << "Number of segregating sites: " << num_segregating << std::endl;
        std::cout << "Number of conserved sites: " << num_same << std::endl;
        std::cout << "Fraction of segragating sites: " << (double)num_segregating/targetsConcat[0].length() << std::endl;
        
        num_segregating = countSegregatingSites(nonTargetConcat);
        num_same = nonTargetConcat[0].length() - num_segregating;
        std::cout << std::endl;
        std::cout << "Alignment not in targets length: " << nonTargetConcat[0].length() << std::endl;
        std::cout << "Number of segregating sites: " << num_segregating << std::endl;
        std::cout << "Number of conserved sites: " << num_same << std::endl;
        std::cout << "Fraction of segragating sites: " << (double)num_segregating/nonTargetConcat[0].length() << std::endl;
        
        
        
    } else {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
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
      
    if (argc != 1) 
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    //opt::utrFile = argv[1];
    //opt::targetLociFile = argv[2];
    
}

// ------------------------- FUNCTIONS -------------------------------------

int countSegregatingSites(const std::vector<string>& utrVector) {
    int num_segregating = 0;
    char c;
    for (string::size_type i = 0; i != utrVector[0].length(); i++) {
        c = utrVector[0][i];
        for (std::vector<string>::size_type j = 1; j != utrVector.size()-1; j++) {
            if (utrVector[j][i] != c) {
                num_segregating++;
                break;
            }    
        }
    }
    return num_segregating;
}


void processAnnotationFile(const string& line, std::map<string,std::vector<string> >& annotationMap) {
    std::vector<string> multipleLoci = split(line, '>'); 
    std::vector<string> locusAndName = split(multipleLoci[1], ':');
    if (locusAndName[1][0] == 'm') {
        std::vector<string> locusScaffoldExons = split(locusAndName[0], '|');
        std::vector<string> locusExons = split(locusScaffoldExons[1], ',');
        // for now, dealing with UTRs composed of multiple exons looks pretty difficult
        // lets ignore them
        if (locusExons.size() == 1) {
            std::vector<string> exonStartEnd = split(locusExons[0], '_');
            std::vector<string> nameVec = split(locusAndName[1], ',');
            std::vector<string> geneNameVec = split(nameVec[0], '.');
            string strand = locusScaffoldExons[2];
            string gene = geneNameVec[0] + ".gene." + geneNameVec[2] + "." + geneNameVec[3] + ".aln"; 
            annotationMap[gene].push_back(locusScaffoldExons[0]); 
            annotationMap[gene].push_back(exonStartEnd[1]);
            annotationMap[gene].push_back(exonStartEnd[2]);
            annotationMap[gene].push_back(strand);
            //locusScaffoldExons[0] + "\t" + exonStartEnd[1] + "\t" + exonStartEnd[2];
        } 
    }    
}

int writeMzebraCoordinatesIntoBedFiles(const std::string& fileRoot, std::map<std::string,std::vector<string> >& annotationMap, const std::vector<int>& targetLoci, std::vector<int>::size_type i, std::ofstream*& targetBed, std::ofstream*& nonTargetBed, int lettersBeforeCutStart, int targetLength, int nextBaseBed) {
        if (annotationMap.count(fileRoot) != 0) {  // If the 3'UTR is confined within a single exon
            string scaffold = annotationMap[fileRoot][0];
            string strand = annotationMap[fileRoot][3];
            if (strand == "+") {
                int utrCutStart = atoi(annotationMap[fileRoot][1].c_str()) + lettersBeforeCutStart;
                int utrStart = atoi(annotationMap[fileRoot][1].c_str());
                int tStart = utrStart + targetLoci[i] - 1;
                int tEnd = tStart + targetLength;
                *targetBed << scaffold << "\t" << tStart << "\t" << tEnd << "\t" << fileRoot << "\t0\t" << strand << std::endl;
                if ((tStart - utrStart + nextBaseBed) > 0) {
                    *nonTargetBed << scaffold << "\t" << utrCutStart + nextBaseBed - 1 << "\t" << tStart << "\t" << fileRoot << "\t0\t" << strand << std::endl;
                }
                nextBaseBed = targetLoci[i] + targetLength - utrCutStart;
            } else if (strand == "-") {
                int utrCutEnd = atoi(annotationMap[fileRoot][2].c_str()) - lettersBeforeCutStart;
                //int utrStart = atoi(annotationMap[fileRoot][1].c_str());
                int utrEnd = atoi(annotationMap[fileRoot][2].c_str());
                int tStart = utrEnd - targetLoci[i] - targetLength;
                int tEnd = utrEnd - targetLoci[i];
                *targetBed << scaffold << "\t" << tStart << "\t" << tEnd << "\t" << fileRoot << "\t0\t" << strand << std::endl;
                *nonTargetBed << scaffold << "\t" << tEnd << "\t" << utrCutEnd - nextBaseBed << "\t" << fileRoot << "\t0\t" << strand << std::endl;
                nextBaseBed = utrCutEnd - tStart;
            }
        } else {
            if (i == 0)
                std::cerr << "Multiple exon 3'UTR: " << fileRoot << " Cannot be included in the BED files for Lake Malawi analysis" << std::endl;
        }
    return nextBaseBed;
}

void writeLastMzebraCoordinateToBed(const std::string& fileRoot, std::map<std::string,std::vector<string> >& annotationMap, const std::vector<int>& targetLoci, std::ofstream*& nonTargetBed, int lettersBeforeCutStart, int nextBaseBed, std::string::size_type l) {
    if (annotationMap.count(fileRoot) != 0) {                   // If the 3'UTR is confined within a single exon
        string scaffold = annotationMap[fileRoot][0];
        string strand = annotationMap[fileRoot][3];
        if (strand == "+") {
            int utrCutStart = atoi(annotationMap[fileRoot][1].c_str()) + lettersBeforeCutStart;
            // int utrStart = atoi(annotationMap[fileRoot][1].c_str());
            // int utrEnd = atoi(annotationMap[fileRoot][2].c_str());
            *nonTargetBed << scaffold << "\t" << utrCutStart + nextBaseBed - 1 << "\t" << utrCutStart + nextBaseBed + l << "\t" << fileRoot << "\t0\t" << strand << std::endl;
        } else if (strand == "-") {
            // int utrStart = atoi(annotationMap[fileRoot][1].c_str());
            int utrCutEnd = atoi(annotationMap[fileRoot][2].c_str()) - lettersBeforeCutStart;
            *nonTargetBed << scaffold << "\t" << utrCutEnd - nextBaseBed - l << "\t" << utrCutEnd - nextBaseBed << "\t" << fileRoot << "\t0\t" << strand << std::endl;
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
