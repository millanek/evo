//
//  process_vcf_stats_testing.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 01/10/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_stats_testing.h"
#include "process_vcf_stats_utils.h"

#define SUBPROGRAM "statsTest"

#define DEBUG 1

static const char *STATS_TEST_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE\n"
"Testing statistics tools; the INPUT_FILE is usually a numerical vector \n"
"\n"
"       -h, --help                              display this help and exit\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "h";

enum { OPT_LDHAT };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string inputFile;
    static string secondFile = "";
    
}

int statsTestingMain(int argc, char** argv) {
    parseStatsTestingOptions(argc, argv);
    // Open connections to read from the vcf and reference genome files
    std::ifstream* inFile = new std::ifstream(opt::inputFile.c_str());
    
    std::vector<double> vc;
    
    std::string line;
    while (getline(*inFile, line)) {
        vc.push_back(convertToDouble(line));
    }
    
    std::cout << "Standard deviation1: " << std_dev(vc) << std::endl;
    
    std::vector<double> vc2;
    if (opt::secondFile != "") {
        std::ifstream* sFile = new std::ifstream(opt::secondFile.c_str());
        while (getline(*sFile, line)) {
            vc2.push_back(convertToDouble(line));
        }
        
        //std::cout << "two sample t-statistic (d = 0): " << two_sample_t(vc, vc2) << std::endl;
        std::cout << "two sample p-value (d != 0): " << two_sample_t(vc, vc2) << std::endl;
        std::cout << "two sample p-value (d > 0.001): " << two_sample_t(vc, vc2, 0.001) << std::endl;
    }
    
    std::cout << "chi_square stat (df=3,chi_sq=3): " << chisq_cdf(3.0, 3.0) << std::endl;
    
    
    std::cout << "Factorial of 5 is: " << factorial(5) << std::endl;
    
    return 0; 
}


void parseStatsTestingOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'h':
                std::cout << STATS_TEST_MESSAGE;
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
        std::cout << "\n" << STATS_TEST_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    
    opt::inputFile = argv[optind++];
    
    if (argc - optind == 1) {
        opt::secondFile = argv[optind++];
    }

}
