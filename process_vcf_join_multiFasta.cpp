//
//  process_vcf_join_multiFasta.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 22/02/2014.
//  Copyright (c) 2014 Milan Malinsky. All rights reserved.
//

#include "process_vcf_join_multiFasta.h"

#define SUBPROGRAM "multi-fasta"

#define DEBUG 1

static const char *MULTI_FASTA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] MULTI_FASTA.fa\n"
"A utility tool for dealing with a multi-fasta file\n"
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


int processMultiFastaMain(int argc, char** argv) {
    parseMultiFastaOptions(argc, argv);
    
    string allSequences = readMultiFastaToOneString(opt::fastaFile, 500000000);
    
    if (opt::out.empty()) {
        print80bpPerLineStdOut(std::cout, allSequences);
    } else {
        string outFN = opt::out + ".joined.fa";
        std::ofstream* outAncestralFile = new std::ofstream(outFN.c_str());
        print80bpPerLineFile(outAncestralFile, allSequences);
    }
    
    return 0;
}

void parseMultiFastaOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::out; break;
            case 'h':
                std::cout << MULTI_FASTA_USAGE_MESSAGE;
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
        std::cout << "\n" << MULTI_FASTA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::fastaFile = argv[optind++];
}