//
//  process_vcf_use_map.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 08/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//


#include "process_vcf_utils.h"
#include "process_vcf_IUPAC.h"
#include "process_vcf_use_map.h"
#include "process_vcf_get_sequences.h"
#include <algorithm>

#define SUBPROGRAM "map"

#define DEBUG 1

static const char *MAP_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] GENOME_SEQUENCE.fa [or VCF_FILE.vcf] LINKAGE_GROUP.order\n"
"Join genomic scaffolds into Linkage Groups, output to STD_OUT\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -v, --vcf-file                          The input is a vcf file rather than a genome file\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


static const char* shortopts = "hv";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "vcf-file", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string inputFile;
    static string LGfile;
    static bool bVCF = false;
}


/*
 TO DO:
 High importance:
 - DEAL WITH scaffold_118 not printing (DONE)
 - Make this work for VCF files (DONE)
 - Make sure all scaffolds are loaded including the last one (DONE)
 - Deal with scaffolds with unknown orientation (DONE)
 - Make sure scaffolds not in any Linkage group are also output (DONE)
 Medium importance:
 - Make sure Linkage groups are output in the right order
 - Look into sorting the VCF in C++ (otherwise igvtools sort seems to work just fine)
 */

int mapMain(int argc, char** argv) {
    parseMapOptions(argc, argv);
    
    // Open connections to read from the vcf and reference genome files
    std::ifstream* LGfile = new std::ifstream(opt::LGfile.c_str());
    
    std::map<string, std::vector<std::vector<string> > > LGmap = loadLinkageGroupMap(LGfile);
    std::cerr << "LG map loaded" << std::endl;
    
    if (opt::bVCF) {
        processVCF(LGmap);
    } else {
        processGenome(LGmap);
    }
    
    return 0;

}

void printProcessedVCFline(std::vector<string>& vcfFields, const bool bInLG, const char thisScaffoldOrientation, const std::string& thisLG, const int& LGsizeUpToHere, const int& thisScaffoldSize, const std::string& line) {
    int pos = atoi(vcfFields[1].c_str());
    vcfFields.erase(vcfFields.begin(), vcfFields.begin()+2);
    if (bInLG) {
        if (thisScaffoldOrientation == '+') {
            std::cout << thisLG << "\t" << LGsizeUpToHere + pos << "\t";
            print_vector_stream(vcfFields, std::cout);
        } else if (thisScaffoldOrientation == '-') {
            std::cout << thisLG << "\t" << LGsizeUpToHere + (thisScaffoldSize - pos + 1) << "\t";
            vcfFields[1] = reverseComplementIUPAC(vcfFields[1]);
            vcfFields[2] = reverseComplementIUPAC(vcfFields[2]);
            print_vector_stream(vcfFields, std::cout);
        }
    } else {
        std::cout << line << std::endl;
    }
}


void processVCF(const std::map<string, std::vector<std::vector<string> > >& LGmap) {
    
    std::cerr << "Processing VCF variant calls from: " << opt::inputFile << std::endl;
    std::cerr << "using linkage groups defined in: " << opt::LGfile << std::endl;
    
    string line;
    std::ifstream* vcfFile = new std::ifstream(opt::inputFile.c_str());
    
    string currentScaffoldNum = "";
    bool bInLG;
    string thisLG;
    int LGsizeUpToHere;
    int thisScaffoldSize;
    char thisScaffoldOrientation;
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#')
            std::cout << line << std::endl;
        else {
            std::vector<string> vcfFields = split(line, '\t');
            if (vcfFields[0] == currentScaffoldNum) {
                printProcessedVCFline(vcfFields, bInLG, thisScaffoldOrientation, thisLG, LGsizeUpToHere, thisScaffoldSize, line);
            } else {
                currentScaffoldNum = vcfFields[0];
                std::vector<string> scaffoldNumVec = split(vcfFields[0], '_');
                bInLG = false;
                thisLG = "";
                for(std::map<string, std::vector<std::vector<string> > >::const_iterator it = LGmap.begin(); it != LGmap.end(); ++it) {
                    // std::cerr << "Building " << it->first << std::endl;
                    LGsizeUpToHere = 0;
                    for (std::vector<std::vector<string> >::size_type i = 0; i != it->second.size(); i++) {
                        // std::cerr << "Adding " << it->second[i][0] << std::endl;
                        if (it->second[i][0] == scaffoldNumVec[1]) {
                            bInLG = true;
                            thisLG = it->first;
                            thisScaffoldSize = atoi(it->second[i][2].c_str());
                            if (it->second[i][1] == "+")
                                thisScaffoldOrientation = '+';
                            else if (it->second[i][1] == "-")
                                thisScaffoldOrientation = '-';
                            break;
                        } else {
                            LGsizeUpToHere += atoi(it->second[i][2].c_str());
                            //std::cerr << "Here: " << it->second[i][0] << std::endl;
                            //std::cerr << "LG:" << it->first << " LGsizeUpToHere: " << LGsizeUpToHere << std::endl;
                        }
                    }
                    if (bInLG) break;
                }
                // std::cerr << "LGsizeUpToHere: " << LGsizeUpToHere << std::endl;
                printProcessedVCFline(vcfFields, bInLG, thisScaffoldOrientation, thisLG, LGsizeUpToHere, thisScaffoldSize, line);
            }
        }
    }
}

void processGenome(const std::map<string, std::vector<std::vector<string> > >& LGmap) {
    
    std::cerr << "Joining genomic scaffolds from: " << opt::inputFile << std::endl;
    std::cerr << "using linkage groups defined in: " << opt::LGfile << std::endl;
    
    string line;
    std::vector<string> currentScaffoldNum;
    string nextScaffoldNum;
    string currentScaffoldString;
    std::map<std::string, std::string> scaffoldStrings;
    std::ifstream* genomeFile = new std::ifstream(opt::inputFile.c_str());
    
    getline(*genomeFile, line);
    currentScaffoldNum = split(line, '_');
    int i = 0;
    while (currentScaffoldString != "") {
        i++;
        currentScaffoldString = readScaffold(genomeFile, nextScaffoldNum);
        if (i % 10 == 0)
            std::cerr << "scaffold " << i << " read" << std::endl;
        /*if (i == 118) {
            std::cerr << "scaffold " << i << " read; length: " << currentScaffoldString.length() << std::endl;
            std::cerr << "currentScaffoldNum[1]: " << currentScaffoldNum[1] << std::endl;
        } */
        scaffoldStrings[currentScaffoldNum[1]] = currentScaffoldString;
        currentScaffoldNum = split(nextScaffoldNum,'_');
    }
    scaffoldStrings[currentScaffoldNum[1]] = currentScaffoldString;
    std::cerr << "scaffold " << currentScaffoldNum[1] << " read" << std::endl;
    
    // std::cerr << "Hello?? " << std::endl;
    std::vector<string> scaffoldsInLinkageGroups;
    string LGString; LGString.reserve(100000000);
    for(std::map<string, std::vector<std::vector<string> > >::const_iterator it = LGmap.begin(); it != LGmap.end(); ++it) {
        LGString = "";
        std::cerr << "Building " << it->first << std::endl;
        for (std::vector<std::vector<string> >::size_type i = 0; i != it->second.size(); i++) {
            std::cerr << "Adding " << it->second[i][0] << " length: " << scaffoldStrings[it->second[i][0]].length() << "(" << it->second[i][2] << ")" << std::endl;
            if (it->second[i][1] == "+") {
                LGString.append(scaffoldStrings[it->second[i][0]]);
            } else if (it->second[i][1] == "-") {
                LGString.append(reverseComplementIUPAC(scaffoldStrings[it->second[i][0]]));
            }
            scaffoldsInLinkageGroups.push_back(it->second[i][0]);
        }
        std::cout << ">" << it->first << std::endl;
        print80bpPerLineStdOut(std::cout, LGString);
    }
    // Print remaining scaffolds (the ones that are not in linkage groups
    for(std::map<std::string, std::string>::iterator it = scaffoldStrings.begin(); it != scaffoldStrings.end(); ++it) {
        if (std::find(scaffoldsInLinkageGroups.begin(), scaffoldsInLinkageGroups.end(), it->first) ==  scaffoldsInLinkageGroups.end()) {
            std::cout << ">scaffold_" << it->first << std::endl;
            print80bpPerLineStdOut(std::cout, it->second);
        }
    }
}

// Load up the annotation file
std::map<string, std::vector<std::vector<string> > > loadLinkageGroupMap(std::ifstream*& LGfile) {
    std::map<string, std::vector<std::vector<string> > > LGmap;
    std::vector<std::vector<string> > scaffoldsInLG;
    string currentLG;
    string line;
    getline(*LGfile, line);
    currentLG = line.substr(1, std::string::npos);
    while (getline(*LGfile, line)) {
        if (line[0] == '>') {
            LGmap[currentLG] = scaffoldsInLG;
            currentLG = line.substr(1, std::string::npos);
            scaffoldsInLG.clear();
        } else {
            std::vector<string> thisScaffold = split(line, '\t');
            scaffoldsInLG.push_back(thisScaffold);
        } 
    }
    return LGmap;
}


void parseMapOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'v': opt::bVCF = true; break;
            case 'h':
                std::cout << MAP_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << MAP_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::inputFile = argv[optind++];
    opt::LGfile = argv[optind++];
}
