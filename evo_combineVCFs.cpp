//
//  evo_combineVCFs.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 27/03/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "evo_combineVCFs.h"
#include "process_vcf_annotation_tools.h"

#define SUBPROGRAM "vcf-comb"

#define DEBUG 1

static const char *VCF_COMB_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE1.vcf VCF_FILE2.vcf REF_1.fa REF2_IN_REF1COORDS.fa MAPPABILITY_UNION.bed\n"
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
    static string mappabilityFile;
    static string out = "Joined_VCF";
}

int VCFcombMain(int argc, char** argv) {
    parseVCFcombOptions(argc, argv);
    const string refAllele = "0/0"; const string altAllele = "1/1";
    string line; // for reading the input files
    
    std::istream* vcfFile1 = createReader(opt::vcfFile1.c_str());
    std::istream* vcfFile2 = createReader(opt::vcfFile2.c_str());
    std::ifstream* refFile1 = new std::ifstream(opt::refFile1.c_str());
    std::ifstream* refFile2 = new std::ifstream(opt::refFile2.c_str());
    
    string extraMaskFileName = opt::out + "_extraMask.bed";
    std::ofstream* extraMaskFile = new std::ofstream(extraMaskFileName.c_str());
    
    
    // Read in the whole reference sequences
    string refSeq1; string refSeq2; refSeq1.reserve(50000000); refSeq2.reserve(50000000);
    getline(*refFile1, line); string chrRef1 = line.substr(1,string::npos);
    getline(*refFile2, line); string chrRef2 = line.substr(1,string::npos);
    assert(chrRef1 == chrRef2);
    
    while (getline(*refFile1, line)) { refSeq1.append(line); }
    while (getline(*refFile2, line)) { refSeq2.append(line); }
    assert(chrRef1.length() == chrRef2.length());
    std::transform(refSeq1.begin(), refSeq1.end(),refSeq1.begin(), ::toupper);
    std::transform(refSeq2.begin(), refSeq2.end(),refSeq2.begin(), ::toupper);
    
    std::map<int,string> VCF1; std::map<int,string> VCF2;
    std::vector<string> header1; std::vector<string> header2;
    std::vector<string> samples1; std::vector<string> samples2;
    std::vector<std::string> fields;
    // Now go through the vcf and add the AA fields
    while (getline(*vcfFile1, line)) {
        if (line[0] == '#' && line[1] == '#') { header1.push_back(line); }
        else if (line[0] == '#' && line[1] == 'C') {
            header1.push_back(line);
            fields = split(line, '\t');
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                samples1.push_back(fields[i]);
            }
        } else {
            fields = split(line, '\t');
            VCF1[atoi(fields[1].c_str())] = line;
        }
    }
    while (getline(*vcfFile2, line)) {
        if (line[0] == '#' && line[1] == '#') { header2.push_back(line); }
        else if (line[0] == '#' && line[1] == 'C') {
            header2.push_back(line); fields = split(line, '\t');
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                samples2.push_back(fields[i]);
            }
        } else {
            fields = split(line, '\t');
            VCF2[atoi(fields[1].c_str())] = line;
        }
    }
    
    std::ifstream* accessibleGenomeBed = new std::ifstream(opt::mappabilityFile.c_str());
    std::cerr << "Loading the accessible genome annotation:" << std::endl;
    AccessibleGenome* ag = new AccessibleGenome(accessibleGenomeBed);
    std::vector<bool> acc(refSeq1.length(),false);
    std::vector<std::vector<int>> thisChrMask = ag->BedFeatureMap.at(chrRef1);
    for (int i = 0; i < thisChrMask[0].size(); i++) {
        for (int j = thisChrMask[0][i]; j < thisChrMask[1][i]; j++) { acc[j] = true; }
    }
    std::cerr << "Done" << std::endl;
    
    
    std::cout << header1[0] << std::endl;
    std::cout << header1[1] << std::endl;
    std::cout << header1[2] << std::endl;
    std::cout << header1.back() << "\t"; print_vector(samples2, std::cout);
    int inMask = 0; int noAlignment = 0; int notVariable = 0; int becomesMultiallelic = 0;
    int vcf1Var = 0; int vcf2Var = 0; int sharedVar = 0; int refDifVar = 0;
    std::cerr << "Length of the chromosome/scaffold: " <<  refSeq1.length() << std::endl;
    
    for (int i = 0; i < refSeq1.length(); i++) {
        int pos = i+1;
       /* if (pos == 7801 || pos == 20205) {
            std::cerr << "VCF1.count("<< pos<<"): " <<  VCF1.count(pos) << "VCF2.count("<< pos<<"): " <<  VCF2.count(pos) << std::endl;
            std::cerr << "refSeq1[i] " <<  refSeq1[i] << "refSeq2[i] " << refSeq2[i] << std::endl;
        }*/
        // if (ag->findIfBPaccessible(chrRef1, i+1) == true) continue; // This assumes a mask file is given
        if (acc[i] == true) { inMask++; continue; }
        if (refSeq1[i] == 'N' || refSeq2[i] == 'N') {
            *extraMaskFile << chrRef1 << "\t" << i << "\t" << i+1 << std::endl;
            noAlignment++; continue;
        }
        if (refSeq1[i] == refSeq2[i]) {
            if (VCF1.count(pos) == 0 && VCF2.count(pos) == 0) { notVariable++; continue; }
            else if (VCF1.count(pos) == 1) {
                if (VCF2.count(pos) == 0) {
                    std::cout << VCF1[pos]; for (int j = 0; j < samples2.size(); j++) { std::cout << "\t" << refAllele; }
                    std::cout << std::endl; vcf1Var++;
                }
                if (VCF2.count(pos) == 1) {
                    fields = split(VCF1[pos], '\t'); char altAllele1 = fields[4][0];
                    fields = split(VCF2[pos], '\t'); char refAllele2 = fields[3][0]; char altAllele2 = fields[4][0];
                    if (refAllele2 != refSeq2[i] && complementIUPAC(refAllele2) == refSeq2[i]) { altAllele2 = complementIUPAC(altAllele2);}
                    if (altAllele1 == altAllele2) {
                        std::vector<string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
                        std::cout << VCF1[pos] << "\t"; print_vector(genotypes, std::cout); sharedVar++;
                    } else {
                        becomesMultiallelic++;
                    } // else this would be multiallelic
                }
            }
            else if (VCF2.count(pos) == 1) {
                fields = split(VCF2[pos], '\t');
                if (fields[3][0] != refSeq2[i] && complementIUPAC(fields[3][0]) == refSeq2[i]) { fields[3][0] = complementIUPAC(fields[3][0]);}
                std::vector<string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
                for (int j = 0; j < NUM_NON_GENOTYPE_COLUMNS; j++) { std::cout << fields[j] << "\t"; }
                for (int j = 0; j < samples1.size(); j++) { std::cout << refAllele << "\t"; }
                print_vector(genotypes,std::cout); vcf2Var++;
            }
        } else // refSeq1[i] != refSeq2[i]
        {
            if (VCF1.count(pos) == 0 && VCF2.count(pos) == 0) {
                std::cout << chrRef1 << "\t" << pos << "\t.\t" << refSeq1[i] << "\t" << refSeq2[i] << "\t1000\tPASS\tAC=" << samples2.size()*2 << "\tGT";
                for (int j = 0; j < samples1.size(); j++) { std::cout << "\t" << refAllele; }
                for (int j = 0; j < samples2.size(); j++) { std::cout << "\t" << altAllele; }
                std::cout << std::endl; refDifVar++;
            }
            else if (VCF1.count(pos) == 1) {
                fields = split(VCF1[pos], '\t'); char altAllele1 = fields[4][0];
                /*if (pos == 7801 || pos == 20205) {
                    std::cerr << "altAllele1 " <<  altAllele1 << "refSeq2[i] " << refSeq2[i] << std::endl;
                }*/
                if (altAllele1 != refSeq2[i]) { becomesMultiallelic++; continue; } // This would be multiallelic
                
                if (VCF2.count(pos) == 0) {
                    std::cout << VCF1[pos]; for (int j = 0; j < samples2.size(); j++) { std::cout << "\t" << refAllele; }
                    std::cout << std::endl; vcf1Var++;
                }
                if (VCF2.count(pos) == 1) {
                    fields = split(VCF2[pos], '\t'); char refAllele2 = fields[3][0]; char altAllele2 = fields[4][0];
                    if (refAllele2 != refSeq2[i] && complementIUPAC(refAllele2) == refSeq2[i]) { altAllele2 = complementIUPAC(altAllele2);}
                    /*if (pos == 7801 || pos == 20205) {
                        std::cerr << "altAllele2 " <<  altAllele1 << "refSeq1[i] " << refSeq2[i] << std::endl;
                    }*/
                    if (altAllele2  == refSeq1[i]) { // The alternative allele in the simDia called VCF needs to match the AstCal allele
                        std::vector<string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
                        std::cout << VCF1[pos] << "\t"; print_vector(genotypes, std::cout); sharedVar++;
                    } else {
                        becomesMultiallelic++;
                    } // else this would be multiallelic
                }
            }
            else if (VCF2.count(pos) == 1) {
                fields = split(VCF2[pos], '\t'); char refAllele2 = fields[3][0]; char altAllele2 = fields[4][0];
                if (refAllele2 != refSeq2[i] && complementIUPAC(refAllele2) == refSeq2[i]) { altAllele2 = complementIUPAC(altAllele2);}
                if (altAllele2  == refSeq1[i]) { // The alternative allele in the simDia called VCF needs to match the AstCal allele
                    std::vector<string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
                    for (int j = 0; j < NUM_NON_GENOTYPE_COLUMNS; j++) { std::cout << fields[j] << "\t"; }
                    for (int j = 0; j < samples1.size(); j++) { std::cout << refAllele << "\t"; }
                    print_vector(genotypes,std::cout); vcf2Var++;
                } else {
                    becomesMultiallelic++;
                }
            }
        }
    }
    std::cerr << "Base categories not resulting in a variant:" << std::endl;
    std::cerr << "inMask\t" << inMask << std::endl;
    std::cerr << "noAlignment\t" << noAlignment << std::endl;
    std::cerr << "notVariable\t" << notVariable << std::endl;
    std::cerr << "becomesMultiallelic\t" << becomesMultiallelic << std::endl;
    std::cerr << "Total non-variant sites:\t" << inMask+noAlignment+notVariable+becomesMultiallelic << std::endl;
    std::cout << std::endl;
    std::cerr << "Base categories resulting in a variant:" << std::endl;
    std::cerr << "Difference between reference sequences:\t" << refDifVar << std::endl;
    std::cerr << "Variant in vcf1:\t" << vcf1Var << std::endl;
    std::cerr << "Variant in vcf2:\t" << vcf2Var << std::endl;
    std::cerr << "Shared vcf1 and vcf2 variant:\t" << sharedVar << std::endl;
    std::cerr << "Total variant sites:\t" << refDifVar+vcf1Var+vcf2Var+sharedVar << std::endl;
    std::cout << std::endl;
    std::cerr << "Total sites:\t" << inMask+noAlignment+notVariable+becomesMultiallelic+refDifVar+vcf1Var+vcf2Var+sharedVar << std::endl;
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
    if (argc - optind < 5) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 5)
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
    opt::refFile2 = argv[optind++];
    opt::mappabilityFile = argv[optind++];
}
