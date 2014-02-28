//
//  process_vcf_mt_sequences.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 30/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_mt_sequences.h"
#include "process_vcf_IUPAC.h"
#include "process_vcf_print_routines.h"

/*
 TO DO:
 High priority:
 - deal with the possibility that a scaffold may not have any variants (DONE - apart from scaffold_0)
 
 */


#define SUBPROGRAM "getMtSeq"

#define DEBUG 1

static const char *GETMTSEQ_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE GENOME_SEQUENCE.fa ANNOTATION.gffExtract\n"
"Obtain full genome sequences from a VCF file (e.g. for multiple alignment and phylogenetic analyses), output to STD_OUT\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the output\n"
"                                               (default: sample ids from the vcf file are used)\n"
"       --LDhat                                 generate output sequences for the LDhat program (compatible with v2.1)\n"
"                                               can be used in conjunction with --split\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hs:";

enum { OPT_LDHAT };

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "LDhat",   no_argument, NULL, OPT_LDHAT },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genomeFile;
    static bool bLDhat = false;
    static bool bWholeGenome = false; // Print the whole genome into one file
    static string sampleNameFile;
    
}



int getMtSeqMain(int argc, char** argv) {
    parseGetMtSeqOptions(argc, argv);
    string vcfFileRoot = stripExtension(opt::vcfFile);    
    
    std::cerr << "Generating the mitochondrial genome sequences using variants from: " << opt::vcfFile << std::endl;
    std::cerr << "and the reference genome: " << opt::genomeFile << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* genomeFile = new std::ifstream(opt::genomeFile.c_str());
    
    // Prepare an object for output files
    string mtFileName = vcfFileRoot + "_mtDNA.fa";
    string mtRefName = vcfFileRoot + "_mtRef.fa";
    std::ofstream* mtFile = new std::ofstream(mtFileName.c_str());    
    
    string line;
    string currentScaffoldNum = "";
    string currentScaffoldReference;
    string::size_type inStrPos = 0;
    std::vector<string> sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);;
    string thisScaffoldName; getline(*genomeFile, thisScaffoldName); thisScaffoldName.erase(0,1);
    unsigned int processedVariantCounter = 0;
    
    std::vector<string::size_type> splits;
    
    
    string sc747Ref; sc747Ref.reserve(20000000);
    string sc2036Ref; sc2036Ref.reserve(20000000);
    string mtRef; mtRef.reserve(20000000);
    while (thisScaffoldName != "") {
        std::cerr << "Starting to read " << thisScaffoldName << std::endl;
        if (thisScaffoldName == "scaffold_747") {
            currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
            thisScaffoldName.erase(0,1);
            sc747Ref.append(currentScaffoldReference);
            mtRef.append(currentScaffoldReference);
        } else if (thisScaffoldName == "scaffold_2036") {
            currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
            thisScaffoldName.erase(0,1);
            sc2036Ref.append(currentScaffoldReference);
            mtRef.append(currentScaffoldReference);
        } else {
            currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
            thisScaffoldName.erase(0,1);
        }
    }
    
    std::vector<string> mtStrings;
    mtStrings.resize(sampleNames.size());
    
    for (std::vector<string>::size_type i = 0; i != mtStrings.size(); i++) {
        mtStrings[i].reserve(20000000);
        mtStrings[i] = "";
    }
    
    string scaffold = "scaffold_747";
    while (getline(*vcfFile, line)) {
        if (line[0] == '#')
            continue;
        else {
            processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            // Also get overall depth for this variant
            std::vector<std::string> info = split(fields[7], ';');
            
            
            if (info[0] != "INDEL" && fields[4].length() == 1) {
                bool printedHet = false;
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    //std::cerr << "Going through genotypes1:" << i << std::endl;
                    //std::cerr << scaffoldStrings.size() << " " << inStrPos << " " << fields[1] << " " << currentScaffoldReference.size() << std::endl;
                    std::vector<string> genotypeFields = split(fields[i], ':');
                    // Unphased vcf
                    std::vector<string> genotype = split(genotypeFields[0], '/');
                    if (opt::bLDhat) {
                        for (int j = 0; j != ((atoi(fields[1].c_str()) - 1)-inStrPos); j++) {
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("0");
                        }
                        if (genotype[0] == "0" && genotype[1] == "0")
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("0");
                        else if (genotype[0] == "1" && genotype[1] == "1")
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("1");
                        else {
                            string ambiguityBase = getAmbiguityCode(fields[3], fields[4]);
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("2");
                        }
                    } else {
                        if (fields[0] == "scaffold_747")
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(mtRef.substr(inStrPos, (atoi(fields[1].c_str()) - 1)-inStrPos));
                        else
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(mtRef.substr(inStrPos, (atoi(fields[1].c_str()) + sc747Ref.length() - 1)-inStrPos));
                        //std::cerr << "Going through genotypes2:" << i << std::endl;
                        if (genotype[0] == "0" && genotype[1] == "0")
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(fields[3]);
                        else if (genotype[0] == "1" && genotype[1] == "1")
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(fields[4]);
                        else {
                            mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(fields[4]);
                            if (!printedHet) {
                                std::cout << "Het in mitochondrial sequence:" << sampleNames[i- NUM_NON_GENOTYPE_COLUMNS] << std::endl;
                                std::cout << line << std::endl;
                                printedHet = true;
                            }
                            //string ambiguityBase = getAmbiguityCode(fields[3], fields[4]);
                            //mtStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(ambiguityBase);
                        }
                    }
                    //std::cerr << "Going through genotypes3:" << i << std::endl;
                }
                if (fields[0] == "scaffold_747")
                    inStrPos = atoi(fields[1].c_str());
                else
                    inStrPos = atoi(fields[1].c_str()) + sc747Ref.length();
                
                
                
#ifdef DEBUG
                if (mtRef[inStrPos-1] != fields[3][0]) {
                    std::cerr << "Error!!! Sequence: " << mtRef[inStrPos-1] << " vcf-ref: " << fields[3][0] << std::endl;
                }
#endif
            }

        }
    }
    
    for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
        mtStrings[i].append(mtRef.substr(inStrPos, string::npos));
    }
#ifdef DEBUG
    if (mtStrings[0].length() != mtRef.length()) {
        std::cerr << "Error!!! Reference scaffold/LG length: " << mtRef.length() << " vcf scaffold length: " << mtStrings[0].length() << std::endl;
    }
#endif
    
    for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
        *mtFile << ">" << sampleNames[i] << std::endl;
        print80bpPerLineFile(mtFile, mtStrings[i]);
    }
    
    
    return 0;
}


void parseGetMtSeqOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 's': arg >> opt::sampleNameFile; break;
            case OPT_LDHAT: opt::bLDhat = true; break;
            case 'h':
                std::cout << GETMTSEQ_USAGE_MESSAGE;
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
        std::cout << "\n" << GETMTSEQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::genomeFile = argv[optind++];
}








 