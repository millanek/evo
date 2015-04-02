//
//  remove_lowercase.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 02/04/2015.
//  Copyright (c) 2015 Milan Malinsky. All rights reserved.
//

#include "remove_lowercase.h"

#define SUBPROGRAM "remove-lowercase"

static const char *REMLC_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] FILTERED_VCF_FILE SAMPLE_SETS.txt\n"
"Remove lowercase characters from a fasta file\n"
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
    static string fastaFile;
    static string out;
}

std::string removeLcLetters(std::ifstream& fastaFile, int bytes) {
    std::string filteredSeq; filteredSeq.reserve(bytes); filteredSeq = "";
    string line;
    while (getline(fastaFile, line)) {
        for (int i = 0; i != line.length(); i++) {
            if (isupper(line[i])) {
                filteredSeq += line[i];
            }
        }
    }
    return filteredSeq;
}



int removeLcMain(int argc, char** argv) {
    parseRemoveLcOptions(argc, argv);
    
    std::ifstream* fastaFile = new std::ifstream(opt::fastaFile.c_str());
    string description; getline(*fastaFile, description);
    
    string finalSequence = removeLcLetters(*fastaFile, 500000000);
    
    if (opt::out.empty()) {
        std::cout << description << std::endl;
        print80bpPerLineStdOut(std::cout, finalSequence);
    } else {
        string outFN = opt::out + ".joined.fa";
        std::ofstream* outFilteredFile = new std::ofstream(outFN.c_str());
        *outFilteredFile << description << std::endl;
        print80bpPerLineFile(outFilteredFile, finalSequence);
    }
    return 0;
}

void parseRemoveLcOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::out; break;
            case 'h':
                std::cout << REMLC_USAGE_MESSAGE;
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
        std::cout << "\n" << REMLC_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::fastaFile = argv[optind++];
}



