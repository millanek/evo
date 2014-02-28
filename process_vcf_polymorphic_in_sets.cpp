//
//  process_vcf_polymorphic_in_sets.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 31/10/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_polymorphic_in_sets.h"
#define SUBPROGRAM "polymorphic"

static const char *POLYMORPHIC_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] FILTERED_VCF_FILE SAMPLE_SETS.txt\n"
"Searching for variants that are polymorphic within particular set(s) of samples\n"
"The SAMPLE_SETS.txt file should have exactly two lines, with each line defining one of the two sample sets\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the VCF file\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "hn:s:";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "run-name",   required_argument, NULL, 'n' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string sampleSets;
    static string sampleNameFile;
    static string runName = "";
}

bool findIfPolymorhicInSet(const std::vector<std::string>& fields, const std::vector<size_t>& set_loci) {
    std::vector<short> setAltAllelePresence;
    short setAltAlleleCount = 0;
    setAltAllelePresence.assign((set_loci.size()),0);
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') {
            if (std::find(set_loci.begin(), set_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set_loci.end()) {
                setAltAllelePresence[i- NUM_NON_GENOTYPE_COLUMNS]++;
                setAltAlleleCount++;
            }
        }
        if (fields[i][2] == '1') {
            if (std::find(set_loci.begin(), set_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set_loci.end()) {
                setAltAllelePresence[i- NUM_NON_GENOTYPE_COLUMNS]++;
                setAltAlleleCount++;
            }
        }
    }
    
    if (setAltAlleleCount == 0) {
        return false;
    } else if (setAltAlleleCount == (set_loci.size() * 2)) {
        return false;
    } else {
        return true;
    }
}

int getNumHets(SetCounts& counts) {
    int num_hets = 0;
    for (std::vector<std::vector<int> >::size_type i = 0; i < counts.individualsWithVariant.size(); i++) {
        if (counts.individualsWithVariant[i] == 1)
            num_hets++;
    }
    return num_hets;
}

int polymorphicMain(int argc, char** argv) {
    parsePolymorphicOptions(argc, argv);
    string fileRoot = stripExtension(opt::sampleSets);
    
    std::cerr << "Filtering a VCF file: " << opt::vcfFile << std::endl;
    std::cerr << "so that only sites that are ploymorphic in sample sets defined in: " << opt::sampleSets << " are output" << std::endl;
    
    // Open connection to read from the vcf file
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::sampleSets.c_str());
    string PolymorphicFileName = fileRoot + "_" + opt::runName + "_polymorphic.vcf";
    std::ofstream* pPolymorphicVCF = new std::ofstream(PolymorphicFileName.c_str());
    
    
    std::vector<std::vector<string> > sets;
    string setString;
    while (getline(*setsFile, setString)) {
        std::vector<string> thisSet = split(setString, ',');
        std::sort(thisSet.begin(),thisSet.end());
        sets.push_back(thisSet);
    }
    
    int numChromosomes;
    int totalVariantNumber = 0;
    string line;
    std::vector<string> sampleNames;
    std::vector<string> fields;
    std::vector<std::vector<size_t> > setsLoci;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            *pPolymorphicVCF << line << std::endl;
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
            
            for (std::vector<std::vector<string> >::size_type i = 0; i != sets.size(); i++) {
                std::vector<size_t> thisSetLoci = locateSet(sampleNames, sets[i]);
                setsLoci.push_back(thisSetLoci);
                std::cerr << "Set" << i << " loci: " << std::endl;
                print_vector_stream(thisSetLoci, std::cerr);
            }
            
            *pPolymorphicVCF << line << std::endl;
        } else {
            totalVariantNumber++;
            
            std::vector<std::string> fields = split(line, '\t');
            
            bool polymorphicInSets = true;
            
            for (std::vector<std::vector<size_t> >::size_type i = 0; i != setsLoci.size(); i++) {
                if(!findIfPolymorhicInSet(fields, setsLoci[i])) {
                    polymorphicInSets = false;
                }
            }
            
            if (polymorphicInSets)
                *pPolymorphicVCF << line << std::endl;
                
        }
    }
    
    
    
    return 0;
}

void parsePolymorphicOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 's': arg >> opt::sampleNameFile; break;
            case 'n': arg >> opt::runName; break;
            case '?': die = true; break;
            case 'h': std::cout << POLYMORPHIC_USAGE_MESSAGE;
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
        std::cout << "\n" << POLYMORPHIC_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::sampleSets = argv[optind++];
}