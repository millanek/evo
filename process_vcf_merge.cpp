//
//  process_vcf_merge.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 20/11/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_merge.h"
#define SUBPROGRAM "merge"

static const char *MERGE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE1.vcf VCF_FILE2.vcf\n"
"Merge two VCF files that contain exactly the same variants (loci), but different samples\n"
"\n"
"       -h, --help                                      display this help and exit\n"
"       -o OUT.vcf, --output-file=OUT.vcf               (required) output file name\n"
"       --genotype-only                                 (optional) include only genotypes from VCF_FILE2.vcf (only the GT field)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1, OPT_GT_ONLY };

static const char* shortopts = "ho:";

static const struct option longopts[] = {
    { "output-file",   required_argument, NULL, 'o' },
    { "samples",   required_argument, NULL, 's' },
    { "genotype-only",   no_argument, NULL, OPT_GT_ONLY },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string VCF1File;
    static string VCF2File;
    static string out;
    static bool bGTonly = false;
}

int mergeMain(int argc, char** argv) {
    parseMergeOptions(argc, argv);
    
    std::cerr << "Merging: " << opt::VCF1File << " and " << opt::VCF2File << std::endl;
    std::cerr << "writing output to:" << opt::out << std::endl;
    
    // Open connection to read from the vcf files
    std::ifstream* pVCF1 = new std::ifstream(opt::VCF1File.c_str());
    std::ifstream* pVCF2 = new std::ifstream(opt::VCF2File.c_str());
    std::ofstream* pOut = new std::ofstream(opt::out.c_str());
    
    
    string line1; string line2;
    while (getline(*pVCF1, line1)) {
        if (line1[1] != '#') break;
        *pOut << line1 << std::endl;
    }
    while (getline(*pVCF2, line2)) {
        if (line2[1] != '#') break;
    }
    
    // Print the sample names
    std::vector<std::string> fields2 = split(line2, '\t');
    std::vector<std::string> genotypes(fields2.begin()+9,fields2.end());
    *pOut << line1 << "\t"; print_vector(genotypes, *pOut);
    
    // Print the variant details
    while (getline(*pVCF1, line1)) {
        getline(*pVCF2, line2);
        fields2 = split(line2, '\t');
        std::vector<std::string> genotypes(fields2.begin()+9,fields2.end());
        if (opt::bGTonly) {
            for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
                std::vector<std::string> gtFields = split(genotypes[i], ':');
                genotypes[i] = gtFields[0];
            }
        }
        *pOut << line1 << "\t"; print_vector(genotypes, *pOut);
    }
    
    return 0;
}

void parseMergeOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'o': arg >> opt::out; break;
            case OPT_GT_ONLY: opt::bGTonly = true; break;
            case '?': die = true; break;
            case 'h': std::cout << MERGE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::out.empty()) {
        std::cerr << "You must specify the name of the output file\n";
        die = true;
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
        std::cout << "\n" << MERGE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::VCF1File = argv[optind++];
    opt::VCF2File = argv[optind++];
}