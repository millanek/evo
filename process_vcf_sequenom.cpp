//
//  process_vcf_sequenom.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 29/08/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_annotation_tools.h"
#include "process_vcf_get_sequences.h"
#include "process_vcf_print_routines.h"
#include "process_vcf_sequenom.h"

#define SUBPROGRAM "sequenom"

#define DEBUG 1

static const char *SEQUENOM_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE GENOME_SEQUENCE.fa ANNOTATION.gffExtract\n"
"Generate a txt file for sequenom assay....\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the output\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hs:";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genomeFile;
    static string geneFile;
    static string sampleNameFile;
    
}

int sequenomMain(int argc, char** argv) {
    parseSequenomOptions(argc, argv);
    std::ios_base::openmode mode = std::ios_base::in;
    string vcfFileName = opt::vcfFile;
    string vcfFileRoot = stripExtension(vcfFileName);
    string genomeFileName = opt::genomeFile;
    string geneFileName = opt::geneFile;
    
    
    std::cerr << "Generating a txt file for sequenom assay from: " << vcfFileName << std::endl;
    std::cerr << "and the reference genome: " << genomeFileName << std::endl;
    std::cerr << "with SNP flanking regions defined in: " << geneFileName << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::ifstream* vcfFile = new std::ifstream(vcfFileName.c_str(), mode);
    std::ifstream* genomeFile = new std::ifstream(genomeFileName.c_str(), mode);
    std::ifstream* geneFile = new std::ifstream(geneFileName.c_str(), mode);
    std::ofstream* sequenomDesignFile = new std::ofstream("sequenom.txt");
    
    string line;
    int currentScaffoldNum = -1;
    string currentScaffoldReference;
    string::size_type inStrPos;
    std::vector<string> sampleNames;
    string thisScaffoldName;
    std::string sequenomString = "";
    
    // Load up the annotation file with SNP flanking regions
    std::vector<std::vector<string> > annotation;
    Annotation wgAnnotation(geneFile, false);
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            
            std::vector<std::string> fields = split(line, '\t');
            // Initialize vectors
            sequenomString.reserve(100000000);
            sequenomString = "";
            
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                if (opt::sampleNameFile.empty())
                    sampleNames.push_back(fields[i]);
                else
                    sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
        } else {
            std::vector<std::string> fields = split(line, '\t');
            // Also get overall depth for this variant
            std::vector<std::string> info = split(fields[7], ';');
            std::vector<std::string> sc = split(fields[0], '_');
            if (atoi(sc[1].c_str()) > currentScaffoldNum) {
                if (currentScaffoldNum >= 0) {
                    // Fill in the rest of the scaffold sequence (after the last variant)
                    sequenomString.append(currentScaffoldReference.substr(inStrPos, string::npos));
                    
#ifdef DEBUG
                    if (sequenomString.length() != currentScaffoldReference.length()) {
                        std::cerr << "Error!!! Reference scaffold length: " << currentScaffoldReference.length() << " vcf scaffold length: " << sequenomString.length() << std::endl;
                    }
#endif
                    
                    string sc = "scaffold_" + numToString(currentScaffoldNum);
                    std::cerr << "Scaffold: " << sc << std::endl;
                    annotation = wgAnnotation.annotationMap[sc]; // Get the SNP surrounding regions on this scaffold
                    for (std::vector<std::vector<string> >::size_type k = 0; k != annotation.size(); k++) {
                        std::vector<string> annotLineVec = split(annotation[k][0], '\t');
                        // std::cerr << "Gene:" << annotLineVec[4] << std::endl;
                        
                        string regionSeq = getIndividualSequenceForThisRegion(annotation[k], annotLineVec[3], sequenomString);
                        string refSeq = getReferenceForThisRegion(annotation[k], annotLineVec[3], currentScaffoldReference);
                        int coordinate = atoi(annotLineVec[2].c_str()) - 100;
                        string location = sc + "_" + numToString(coordinate);
                        
                        writeSequenomOutput(regionSeq, refSeq, location, sequenomDesignFile);
                    }
                    
                    sequenomString = "";
                } else {
                    getline(*genomeFile, thisScaffoldName);
                }
                
                std::cerr << "Starting to read " << thisScaffoldName << std::endl;
                currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
#ifdef DEBUG2
                std::cerr << "Finished reading. Next scaffold will be: " << thisScaffoldName << std::endl;
#endif
                inStrPos = 0;
                currentScaffoldNum = atoi(sc[1].c_str());
            }
            if (info[0] != "INDEL") {
                sequenomString.append(currentScaffoldReference.substr(inStrPos, (atoi(fields[1].c_str()) - 1)-inStrPos));
                Counts thisVariantCounts = getThisVariantCounts(fields);
                if (thisVariantCounts.overall != ((fields.size()-NUM_NON_GENOTYPE_COLUMNS)*2)) {
                    string ambiguityBase = getAmbiguityCode(fields[3], fields[4]);
                    sequenomString.append(ambiguityBase);
                } else {
                    sequenomString.append(fields[4]);
                }
                inStrPos = atoi(fields[1].c_str());
                
#ifdef DEBUG
                if (currentScaffoldReference[inStrPos-1] != fields[3][0]) {
                    std::cerr << "Error!!! Sequence: " << currentScaffoldReference[inStrPos-1] << " vcf-ref: " << fields[3][0] << std::endl;
                }
#endif
            }
        }
    }
    
    return 0;
}


bool checkAllDNA(const std::string& seq) {
    //std::cerr << seq << std::endl;
    bool bDNA = true;
    for (std::string::size_type i = 0; i != seq.length(); i++) {
        if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
            bDNA = false;
        }
    }
    return bDNA;
}

void writeSequenomOutput(const std::string& regionSeq, const std::string& refSeq, const std::string& loc, std::ofstream*& sequenomFile) {
    
    char refBase = refSeq[100];
    char altBase = regionSeq[100];
    switch (altBase)
    {
        case 'K': altBase = (refBase == 'G') ? 'T' : 'G'; break;
        case 'M': altBase = (refBase == 'A') ? 'C' : 'A'; break;
        case 'R': altBase = (refBase == 'A') ? 'G' : 'A'; break;
        case 'S': altBase = (refBase == 'C') ? 'G' : 'C'; break;
        case 'W': altBase = (refBase == 'A') ? 'T' : 'A'; break;
        case 'Y': altBase = (refBase == 'C') ? 'T' : 'C'; break;
        default: altBase = 'X'; refBase = 'X'; break;
    }
    
    bool bPrimerCanBeDesigned = true;
    bool bRightAllDNA = checkAllDNA(regionSeq.substr(101,28));
    bool bLeftAllDNA = checkAllDNA(regionSeq.substr(72,28));
    if (bRightAllDNA == false && bLeftAllDNA == false) {
        bPrimerCanBeDesigned = false;
    }
    
    if (bPrimerCanBeDesigned)
        *sequenomFile << loc << "\t" << regionSeq.substr(0,100) << "[" << refBase << "/" << altBase << "]" << regionSeq.substr(101, string::npos) << std::endl;
    else
        *sequenomFile << loc << "\t" << "Primer problem 28bp" << std::endl;
    
}

void parseSequenomOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 's': arg >> opt::sampleNameFile; break;
            case 'h':
                std::cout << SEQUENOM_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 3) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << SEQUENOM_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::genomeFile = argv[optind++];
    opt::geneFile = argv[optind++];
}


