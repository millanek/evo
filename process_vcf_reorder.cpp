//
//  process_vcf_reorder.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 23/10/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_reorder.h"
#define SUBPROGRAM "reorder"

static const char *REORDER_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE.vcf NEW_ORDER_FILE.txt\n"
"\n"
"The SAMPLE_SETS.txt file should have exactly two lines, with each line defining one of the two sample sets\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the VCF file\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "hs:n:";

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "samples",   required_argument, NULL, 's' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string newOrderFile;
    static string sampleNameFile;
    static string runName;
}

std::map<string,size_t> linkVectors(const std::vector<string>& newOrder, const std::vector<string> oldOrder) {
    std::map<string,size_t> link;
    for (std::vector<string>::size_type i = 0; i != newOrder.size(); i++) {
        std::vector<string>::const_iterator it = std::find(oldOrder.begin(), oldOrder.end(), newOrder[i]);
        if (it == oldOrder.end()) {
            if (opt::sampleNameFile.empty())
                std::cerr << "The column names in NEW_ORDER_FILE.txt do not correspond to names in the VCF header" << std::endl;
            else
                std::cerr << "The column names in NEW_ORDER_FILE.txt do not correspond to names in " << opt::sampleNameFile << std::endl;
            exit(1);
        }
        size_t pos = std::distance(oldOrder.begin(), it);
        link[oldOrder[i]] = pos;
    }
    return link;
}


int reorderMain(int argc, char** argv) {
    parseReorderOptions(argc, argv);
    string fileRoot = stripExtension(opt::vcfFile);
    
    std::cerr << "Reordering columns in: " << opt::vcfFile << std::endl;
    std::cerr << "using ordering in: " << opt::newOrderFile << std::endl;
    
    // Open connection to read from the vcf file
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    string reorderedFileName = fileRoot + opt::runName + "_reordered.vcf";
    std::ofstream* pReordered = new std::ofstream(reorderedFileName.c_str());
    
    int numChromosomes;
    int totalVariantNumber = 0;
    string line;
    std::vector<string> sampleNames;
    std::vector<string> newOrder = readSampleNamesFromTextFile(opt::newOrderFile);
    std::vector<string> fields;
    std::map<string, size_t> link;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            *pReordered << line << std::endl;
        } else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
            numChromosomes = (int)numSamples * 2;
            // std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            assert(sampleNames.size() == newOrder.size());
            link = linkVectors(newOrder, sampleNames);
            for (std::vector<std::string>::size_type i = 0; i != NUM_NON_GENOTYPE_COLUMNS; i++) {
                *pReordered << fields[i] << "\t";
            }
            
            for (std::vector<std::string>::size_type i = 0; i != sampleNames.size() - 1; i++) {
                *pReordered << newOrder[i] << "\t";
            } *pReordered << newOrder[newOrder.size()-1] << std::endl;
        } else {
            totalVariantNumber++;
            
            std::vector<std::string> fields = split(line, '\t');
            for (std::vector<std::string>::size_type i = 0; i != NUM_NON_GENOTYPE_COLUMNS; i++) {
                *pReordered << fields[i] << "\t";
            }
            
            for (std::vector<std::string>::size_type i = 0; i != sampleNames.size() - 1; i++) {
                *pReordered << fields[link[sampleNames[i]]+NUM_NON_GENOTYPE_COLUMNS] << "\t";
            } *pReordered << fields[link[sampleNames[sampleNames.size()-1]]+NUM_NON_GENOTYPE_COLUMNS] << std::endl;
 
            
        }
    }
    
    
    
    return 0;
}

void parseReorderOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'n': arg >> opt::runName; break;
            case 's': arg >> opt::sampleNameFile; break;
            case '?': die = true; break;
            case 'h': std::cout << REORDER_USAGE_MESSAGE;
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
        std::cout << "\n" << REORDER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::newOrderFile = argv[optind++];
}
