//
//  evo_combineVCFs.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 27/03/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "evo_combineVCFs.h"

#define SUBPROGRAM "vcf-comb"

#define DEBUG 1

static const char *VCF_COMB_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE1.vcf VCF_FILE2.vcf REF_1.fa REF2_IN_REF1COORDS.fa\n"
"Generates a joint VCF file from two VCFs that have different samples and were originally called against two different references\n"
"This is done per-chromosome\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -o, --out=FILE_ROOT                     the output file will be 'FILE_ROOT_Joint.vcf'\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


static const char* shortopts = "ho:";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "out",   required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile1;
    static string vcfFile2;
    static string refFile1;
    static string refFile2;
    static string out = "";
}

int VCFcombMain(int argc, char** argv) {
    parseVCFcombOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* vcfFile1 = createReader(opt::vcfFile1.c_str());
    std::istream* vcfFile2 = createReader(opt::vcfFile2.c_str());
    std::ifstream* refFile1 = new std::ifstream(opt::refFile1.c_str());
    std::ifstream* refFile2 = new std::ifstream(opt::refFile2.c_str());
    
    // Read in the whole reference sequences
    string refSeq1; string refSeq2; refSeq1.reserve(50000000); refSeq2.reserve(50000000);
    getline(*refFile1, line); string chrRef1 = line.substr(1,string::npos);
    getline(*refFile2, line); string chrRef2 = line.substr(1,string::npos);
    assert(chrRef1 == chrRef2);
    
    while (getline(*refFile1, line)) { refSeq1.append(line); }
    while (getline(*refFile2, line)) { refSeq2.append(line); }
    assert(chrRef1.length() == chrRef2.length());
    
    std::map<int,string> VCF1; std::map<int,string> VCF2;
    std::vector<string> header1; std::vector<string> header2;
    std::vector<string> samples1; std::vector<string> samples2;
    std::vector<std::string> fields;
    // Now go through the vcf and add the AA fields
    while (getline(*vcfFile1, line)) {
        if (line[0] == '#' && line[1] == '#') { header1.push_back(line); }
        else if (line[0] == '#' && line[1] == 'C') { samples1 = split(line, '\t'); }
        else {
            fields = split(line, '\t');
            VCF1[atoi(fields[1].c_str())] = line;
        }
    }
    while (getline(*vcfFile2, line)) {
        if (line[0] == '#' && line[1] == '#') { header2.push_back(line); }
        else if (line[0] == '#' && line[1] == 'C') { samples2 = split(line, '\t'); }
        else {
            fields = split(line, '\t');
            VCF2[atoi(fields[1].c_str())] = line;
        }
    }
    
    return 0;
}


void parseVCFcombOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::out; break;
            case 'h':
                std::cout << VCF_COMB_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 4) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 4)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << VCF_COMB_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile1 = argv[optind++];
    opt::vcfFile2 = argv[optind++];
    opt::refFile1 = argv[optind++];
    opt::refFile1 = argv[optind++];
}
