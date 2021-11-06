//
//  evo_agpToNewFasta.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 03/11/2021.
//  Copyright © 2021 Milan Malinsky. All rights reserved.
//

#include "process_vcf_utils.h"
#include "evo_agpToNewFasta.h"



#define SUBPROGRAM "agpToNewFasta"

#define DEBUG 1
#define HUNDRED_Ns "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

static const char *AGPTOFASTA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] AGP_FILE.agp GENOME_SEQUENCE.fa\n"
"Obtain a new genome sequence with rearrangements as indicated in an AGP file; output to STD_OUT\n"
"\n"
"       -h, --help                                  display this help and exit\n"
"       --by-scaffold                               output by scaffold/LG (each scaffold/LG) has its own file with sequences\n"
"                                                   for all samples\n"
"       --whole-genome                              output is one file with the whole genome concatenated for all samples\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_BY_SCAFFOLD, OPT_WG };

static const char* shortopts = "hs:H:";

static const struct option longopts[] = {
    { "by-scaffold",   no_argument, NULL, OPT_BY_SCAFFOLD },
    { "whole-genome",   no_argument, NULL, OPT_WG },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string agpFile;
    static string genomeFile = "";
    static string outgroupFile;
    static bool bByScaffold = false;
    static bool bWholeGenome = false; // Print the whole genome into one file

}


std::map<string, std::vector<std::vector<string> > > readAGPstructureToMap(std::istream* agpFile) {
    std::map<string, std::vector<std::vector<string> > > newAgpStructure;
    string line;
    while (getline(*agpFile, line)) {
        if (line[0] == '#') { continue; } // Skip over any comment lines
        else {
            std::vector<std::string> fields = split(line, '\t');
            string currentScaffold = fields[0];
            newAgpStructure[currentScaffold].push_back(fields);
        }
    }
    return newAgpStructure;
}

int agpFastaMain(int argc, char** argv) {
    
    parseAgpFastaOptions(argc, argv);
    
    std::cerr << "Generating a new genome sequence with rearrangements as indicated in an AGP file: " << opt::agpFile << std::endl;
    std::cerr << "from the original reference genome: " << opt::genomeFile << std::endl;
    
    // Read from the agp and reference genome files
    std::istream* agpFile = createReader(opt::agpFile.c_str());
    std::map<string, string> fastaSeqs = readMultiFastaToMap(opt::genomeFile);
    std::map<string, std::vector<std::vector<string> > > newAgpStructure = readAGPstructureToMap(agpFile);
    std::map<string, string> outputSeqs;
    
    for(std::map<string,std::vector<std::vector<string> > >::iterator it = newAgpStructure.begin(); it != newAgpStructure.end(); ++it) {
        string processedScaffold = it->first;
        std::vector<std::vector<string> > processedScaffoldStructure = it->second;
        outputSeqs[processedScaffold] = ""; outputSeqs[processedScaffold].reserve(50000000);
        
        for (std::vector<std::string>::size_type i = 0; i != processedScaffoldStructure.size(); i++) {
            string toAdd;
            if (processedScaffoldStructure[i][6] == "scaffold") {
                toAdd = HUNDRED_Ns;
            } else {
                string originScaffold = processedScaffoldStructure[i][5];
                int originStart = atoi(processedScaffoldStructure[i][6].c_str()); int originEnd = atoi(processedScaffoldStructure[i][7].c_str());
                int originLength = originEnd - originStart;
                toAdd = fastaSeqs.at(originScaffold).substr(originStart-1,originLength);
                // If the sequence should be inverted, get its reverse complement
                if (processedScaffoldStructure[i][8] == "-") {
                    std::reverse(toAdd.begin(), toAdd.end());
                    // complement
                    for (std::string::size_type j = 0; j != toAdd.length(); j++) {
                        switch (toAdd[j]) {
                            case 'A': toAdd[j] = 'T'; break;
                            case 'T': toAdd[j] = 'A'; break;
                            case 'C': toAdd[j] = 'G'; break;
                            case 'G': toAdd[j] = 'C'; break;
                            default:
                                break;
                        }
                    }
                }
            }
            outputSeqs[processedScaffold].append(toAdd);
        }
    }
    
    for(std::map<string,std::vector<std::vector<string> > >::iterator it = newAgpStructure.begin(); it != newAgpStructure.end(); ++it) {
        string processedScaffold = it->first;
        std::cerr << "New scaffold " << processedScaffold << " length: " << outputSeqs[processedScaffold].length() << std::endl;
        std::cout << ">" << processedScaffold << std::endl;
        print80bpPerLineStdOut(std::cout, outputSeqs[processedScaffold]);
    }
    
    return 0;
}



void parseAgpFastaOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'h':
                std::cout << AGPTOFASTA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    
    if (die) {
        std::cout << "\n" << AGPTOFASTA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::agpFile = argv[optind++];
    opt::genomeFile = argv[optind++];
}

