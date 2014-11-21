//
//  process_vcf_shortRNA.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 03/03/2014.
//  Copyright (c) 2014 Milan Malinsky. All rights reserved.
//

#include "process_vcf_shortRNA.h"

#define SUBPROGRAM "smallRNA"

#define DEBUG 1

static const char *SMALL_RNA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] shortRNAfile.fa\n"
"Generate a distribution showing read lengths and the starting nucleotide for a smallRNA library\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -o, --out=FILE_ROOT                     the outpur file will be 'FILE_ROOT.ancestralSequence.fa'\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "ho:";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "out",   no_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string smallRNAFile;
    static string out;
    static bool bAncFromMaf = false;
}


int shortRNAMain(int argc, char** argv) {
    parseShortRNAOptions(argc, argv);
    
    std::ifstream* smallRNAFile = new std::ifstream(opt::smallRNAFile);
    
    string refFastaFileRoot;
    if (opt::out.empty()) {
        refFastaFileRoot = stripExtension(opt::smallRNAFile);
    } else {
        refFastaFileRoot = opt::out;
    }
    
    string outDistribution = refFastaFileRoot + "smallRNAdist.forR";
    string outDistributionUnique = refFastaFileRoot + "smallRNAdistUnique.forR";
  
    std::cerr << "Output will be in: " << outDistribution << std::endl;
    std::cerr << "and: " << outDistributionUnique << std::endl;
    std::ofstream* outDistributionFile = new std::ofstream(outDistribution);
    std::ofstream* outDistributionUniqueFile = new std::ofstream(outDistributionUnique);
    
    std::map<string, bool> seen;
    
    string line;
    getline(*smallRNAFile, line);
    bool fastq = false;
    if (line[0] == '@') { fastq = true; } else if (line[0] == '>') { fastq = false; } else { std::cerr << "File format not recognised" << std::endl; exit(1); }
    int linesRead = 1;
    
    std::vector<std::vector<int> > counts; std::vector<std::vector<int> > uniqueCounts; counts.resize(34); uniqueCounts.resize(34);
    for (std::vector<std::vector<int> >::size_type i = 0; i != counts.size(); i++) {
        counts[i].resize(4); uniqueCounts[i].resize(4);
    }
    
    while (getline(*smallRNAFile, line)) {
        linesRead++;
        if ((!fastq && (linesRead % 2 == 0))
            || (fastq && ((linesRead+2) % 4 == 0))) {
            if (line[0] == 'N') { continue; }
            else {
                // std::cerr << "Seq: " << line << std::endl;
                string::size_type l = line.length();
                if (seen.find(line) == seen.end()) {
                    switch (line[0]) {
                        case 'A': counts[l][0]++; uniqueCounts[l][0]++;
                            // std::cerr << "l: " << l << " fb: " << line[0] << "c: " << counts[l][0] << std::endl;
                            break;
                        case 'C': counts[l][1]++; uniqueCounts[l][1]++; break;
                        case 'G': counts[l][2]++; uniqueCounts[l][2]++; break;
                        case 'T': counts[l][3]++; uniqueCounts[l][3]++; break;
                        default: std::cerr << "Unknown base in the sequence: " << line[0] << " Exiting..." << std::endl; exit(1);
                    }
                    seen[line] = true;
                } else {
                    switch (line[0]) {
                        case 'A': counts[l][0]++;
                            // std::cerr << "l: " << l << " fb: " << line[0] << "c: " << counts[l][0] << std::endl;
                            break;
                        case 'C': counts[l][1]++; break;
                        case 'G': counts[l][2]++; break;
                        case 'T': counts[l][3]++; break;
                        default: std::cerr << "Unknown base in the sequence: " << line[0] << " Exiting..." << std::endl; exit(1);
                    }
                }
            }
        }
    }
    
    //std::cerr << "Counts_size: " << counts.size() << " " << counts[18][1] << std::endl;
    
    for (std::vector<int>::size_type j = 0; j != 4; j++) {
        for (std::vector<std::vector<int> >::size_type i = 18; i != counts.size(); i++) {
            //std::cerr << " i: " << i << " " << counts[i][j] << std::endl;
            if (i == (counts.size()-1)) {
                *outDistributionFile << counts[i][j] << std::endl;
                *outDistributionUniqueFile << uniqueCounts[i][j] << std::endl;
            }
            else {
                *outDistributionFile << counts[i][j] << "\t";
                *outDistributionUniqueFile << uniqueCounts[i][j] << "\t";
            }
        }
    }
    
    return 0;
}

void parseShortRNAOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::out; break;
            case 'h':
                std::cout << SMALL_RNA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << SMALL_RNA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::smallRNAFile = argv[optind++];
}



