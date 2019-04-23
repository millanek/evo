//
//  process_vcf_fill_aa.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 18/02/2014.
//  Copyright (c) 2014 Milan Malinsky. All rights reserved.
//

#include "process_vcf_fill_aa.h"

#define SUBPROGRAM "aa-fill"

#define DEBUG 1

static const char *FILL_AA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE.vcf ANCESTRAL_SEQ_IN_REF_COORDS.fa\n"
"Adds the AA (ancestral allele) info field to a VCF file given an ancestral sequence in ref genome coordinates\n"
"for a scaffold\n"
"ANCESTRAL_SEQ_IN_REF_COORDS.fa can be produced by " PROGRAM_BIN " aa-seq\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -o, --out=FILE_ROOT                     the output file will be 'FILE_ROOT_AAfilled.vcf'\n"
"       -i, --addAsAnIndividual=IndividualName  Instead of filling the ancestral allele, this adds the sequence as"
"                                               an individual to the VCF file\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


static const char* shortopts = "ho:i:";

static const struct option longopts[] = {
    { "addAsAnIndividual",   required_argument, NULL, 'i' },
    { "help",   no_argument, NULL, 'h' },
    { "out",   required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string ancSeqFile;
    static string out = "";
    static string IndividualName = "";
}


int fillAaMain(int argc, char** argv) {
    parseFillAaOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* ancSeqFile = new std::ifstream(opt::ancSeqFile.c_str());
    
    string refFastaFileRoot;
    if (opt::out.empty()) {
        refFastaFileRoot = stripExtension(opt::vcfFile);
    } else {
        refFastaFileRoot = opt::out;
    }
    // the output file
    std::ostream* outFile;
    if (opt::out == "") {
        outFile = &std::cout;
    } else {
        string outFN = refFastaFileRoot + "_AAfilled.vcf";
        outFile = createWriter(outFN.c_str());
    }
    
    // Read in the whole ancestral sequence
    std::map<string, string> ancSeqs;
    getline(*ancSeqFile, line);
    string currentScaffold = line.substr(1,string::npos);
    ancSeqs[currentScaffold] = ""; ancSeqs[currentScaffold].reserve(50000000);
    while (getline(*ancSeqFile, line)) {
        if (line[0] != '>') {
            ancSeqs[currentScaffold].append(line);
        } else {
            // std::cerr << currentScaffold << " length: " << ancSeqs[currentScaffold].length() << std::endl;
            currentScaffold = line.substr(1,string::npos);
            ancSeqs[currentScaffold] = ""; ancSeqs[currentScaffold].reserve(50000000);
        }
    }
    
    // Now go through the vcf and add the AA fields
    int totalVariantNumber = 0;
    int aaDashCount = 0; int aaRefCount = 0; int aaAltCount = 0; int aaDiffCount = 0; int aaNcount = 0;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            *outFile << line << std::endl;
        else if (line[0] == '#' && line[1] == 'C') {
            if (opt::IndividualName == "") {
                *outFile << "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">" << std::endl;
                *outFile << line << std::endl;
            } else {
                *outFile << line << "\t" << opt::IndividualName << std::endl;
            }
        } else {
            totalVariantNumber++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (info[0] != "INDEL") {
                assert(ancSeqs.find(fields[0]) != ancSeqs.end());
                char AA;
                if (ancSeqs[fields[0]].length() == 0) {
                    AA = 'N';
                } else {
                    AA = ancSeqs[fields[0]][atoi(fields[1].c_str())-1];
                    if (AA == '-') { aaDashCount++; }
                    else if (AA == 'N') { aaNcount++; }
                    else if (AA == fields[3][0]) { aaRefCount++; }
                    else if (AA == fields[4][0]) { aaAltCount++; }
                    else if (!((AA == fields[3][0]) || (AA == fields[4][0]))) {
                        aaDiffCount++;
                        // std::cerr << fields[0] << "\t" << fields[1] << "\t" << fields[3] << "\t" << fields[4] << "\t" << AA << std::endl;
                    }
                    // assert((AA == fields[3][0]) || (AA == fields[4][0]));
                }
                if (opt::IndividualName == "") {
                    fields[7] += ";AA="; fields[7] += AA;
                    print_vector(fields, *outFile, '\t');
                } else {
                    string genotypeToAdd;
                    if (AA == fields[3][0])
                        genotypeToAdd = "0/0";
                    else if (AA == fields[4][0])
                        genotypeToAdd = "1/1";
                    else
                        genotypeToAdd = "./.";
                    fields.push_back(genotypeToAdd);
                    print_vector(fields, *outFile, '\t');
                }
            } else {
                if (opt::IndividualName == "") {
                    *outFile << line << std::endl;
                } else {
                    fields.push_back("./.");
                    print_vector(fields, *outFile, '\t');
                }
            }
            if (totalVariantNumber % 100000 == 0) {
                double totalAAfilled = aaRefCount + aaAltCount + aaDashCount + aaDiffCount + aaNcount;
                std::cerr << totalVariantNumber << " variants processed. AA=Ref:" << aaRefCount << "("<< 100*(aaRefCount/totalAAfilled) <<"%); AA=Alt:" << aaAltCount << "("<< 100*(aaAltCount/totalAAfilled) <<"%); AA='-':" << aaDashCount << "("<< 100*(aaDashCount/totalAAfilled) << "%); AA=?(Neither Ref nor Alt):" << aaDiffCount << "("<< 100*(aaDiffCount/totalAAfilled) <<"%); AA=N:" << aaNcount << "("<< 100*(aaNcount/totalAAfilled) << "%)" << std::endl;
            }
        }
    }
    // Final summary
    double totalAAfilled = aaRefCount + aaAltCount + aaDashCount + aaDiffCount;
    std::cerr << std::endl;
    std::cerr << "All " << totalVariantNumber << " variants processed. AA=Ref:" << aaRefCount << "("<< 100*(aaRefCount/totalAAfilled) <<"%); AA=Alt:" << aaAltCount << "("<< 100*(aaAltCount/totalAAfilled) <<"%); AA='-':" << aaDashCount << "("<< 100*(aaDashCount/totalAAfilled) << "%); AA=?(Neither Ref nor Alt):" << aaDiffCount << "("<< 100*(aaDiffCount/totalAAfilled) <<"%)" << std::endl;
    
    return 0;
}

void parseFillAaOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::out; break;
            case 'i': arg >> opt::IndividualName; break;
            case 'h':
                std::cout << FILL_AA_USAGE_MESSAGE;
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
        std::cout << "\n" << FILL_AA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::ancSeqFile = argv[optind++];
}
