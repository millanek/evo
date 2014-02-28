//
//  process_vcf_abba_baba.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 12/12/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_abba_baba.h"
#include "process_vcf_stats_utils.h"


#define SUBPROGRAM "abba-baba"

#define DEBUG 1

static const char *ABBA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf SETS.txt\n"
"Calculate the D-statistic (abba/baba) as definded in Durand et al. 2011"
"This is a Four-Taxon Statistic to Test for Admixture using data from a VCF file\n"
"The SETS.txt should have exactly four lines:\n"
"Line 1: Outgroup individuals\n"
"Line 2,3,4: P3,P2,P1 individuals respectively (as defined in the Durand et al. 2011 paper\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -f, --frequency                         use allele frequency data instead of single sequences for each of (P1,P2,P3,O)\n"
"       --AAeqO                                 ancestral allele infor in the VCF is from the outgroup (e.g. Pnyererei for Malawi)\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   (optional) supply a file of sample identifiers\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_AA_EQ_O };

static const char* shortopts = "hs:f";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "AAeqO",   no_argument, NULL, OPT_AA_EQ_O },
    { "frequency",   no_argument, NULL, 'f' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string sampleNameFile;
    static bool bFrequency = false;
    static bool bAaEqO = false;
    static int minScLength = 0;
}


void singleSequence() {
    string line; // for reading the input files
    
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    
    //
    string outgroup; size_t Opos; if (!opt::bAaEqO) { getline(*setsFile, outgroup); } else { outgroup = "VCF AA field"; }
    string P3; getline(*setsFile, P3); size_t P3pos;
    string P2; getline(*setsFile, P2); size_t P2pos;
    string P1; getline(*setsFile, P1); size_t P1pos;
    
    // Now go through the vcf and find the ABBA and BABA sites
    int totalVariantNumber = 0;
    int ABBA = 0; int BABA = 0; int XXBA = 0;
    int Dnumerator = 0; int Ddenominator = 0;
    std::vector<string> sampleNames;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            if (opt::bAaEqO) { Opos = locateOneSample(sampleNames, outgroup); } else { Opos = 1000; }
            P3pos = locateOneSample(sampleNames, P3);
            P2pos = locateOneSample(sampleNames, P2); P1pos = locateOneSample(sampleNames, P1);
            
            
            if (! opt::bAaEqO) { assert(fields[Opos+NUM_NON_GENOTYPE_COLUMNS] == outgroup); std::cerr << "Outgroup: " << outgroup << " Pos: " << Opos << std::endl; }
            else { std::cerr << "Outgroup: " << outgroup << std::endl; }
            assert(fields[P3pos+NUM_NON_GENOTYPE_COLUMNS] == P3); std::cerr << "P3: " << P3 << " Pos: " << P3pos << std::endl;
            assert(fields[P2pos+NUM_NON_GENOTYPE_COLUMNS] == P2); std::cerr << "P2: " << P2 << " Pos: " << P2pos << std::endl;
            assert(fields[P1pos+NUM_NON_GENOTYPE_COLUMNS] == P1); std::cerr << "P1: " << P1 << " Pos: " << P1pos << std::endl;
        } else {
            totalVariantNumber++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (info[0] != "INDEL") {
                string AA = split(info[info.size()-1],'=')[1];
                if ((AA == fields[3]) &&
                    (fields[Opos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "0/0") && (fields[P3pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "1/1")) {
                    XXBA++;
                    //
                    if ((fields[P2pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "1/1") && (fields[P1pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "0/0")) {
                        ABBA++; Ddenominator++; Dnumerator++;
                    }
                    // BABA
                    if ((fields[P2pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "0/0") && (fields[P1pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "1/1")) {
                        BABA++; Ddenominator++; Dnumerator--;
                    }
                }
                if ((AA == fields[4]) &&
                    (fields[Opos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "1/1") && (fields[P3pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "0/0")) {
                    XXBA++;
                    // ABBA
                    if ((fields[P2pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "0/0") && (fields[P1pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "1/1")) {
                        ABBA++; Ddenominator++; Dnumerator++;
                    }
                    // BABA
                    if ((fields[P2pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "1/1") && (fields[P1pos+NUM_NON_GENOTYPE_COLUMNS].substr(0,3) == "0/0")) {
                        BABA++; Ddenominator++; Dnumerator--;
                    }
                }
            } else {
                // Ignoring INDELS for now
            }
            if (totalVariantNumber % 100000 == 0) {
                std::cerr << totalVariantNumber << " variants processed. XXBA=" << XXBA << "; ABBA=" << ABBA << "; BABA=" << BABA << "; D=" << (double)Dnumerator/Ddenominator << std::endl;
            }
        }
    }
    // Final summary
    // double totalAAfilled = aaRefCount + aaAltCount + aaDashCount + aaDiffCount;
    std::cerr << std::endl;
    std::cerr << totalVariantNumber << " variants processed. ABBA=" << ABBA << "; BABA=" << BABA << "; D=" << (double)Dnumerator/Ddenominator << std::endl;
    // std::cerr << "All " << totalVariantNumber << " variants processed. AA=Ref:" << aaRefCount << "("<< 100*(aaRefCount/totalAAfilled) <<"%); AA=Alt:" << aaAltCount << "("<< 100*(aaAltCount/totalAAfilled) <<"%); AA='-':" << aaDashCount << "("<< 100*(aaDashCount/totalAAfilled) << "%); AA=?(Neither Ref nor Alt):" << aaDiffCount << "("<< 100*(aaDiffCount/totalAAfilled) <<"%)" << std::endl;
}


void frequencyData() {
    string line; // for reading the input files
    
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    
    //
    string outgroupString; std::vector<size_t> Opos; std::vector<string> outgroup;
    if (!opt::bAaEqO) { getline(*setsFile, outgroupString); outgroup = split(outgroupString, ','); } else { outgroupString = "VCF AA field"; }
    string P3string; getline(*setsFile, P3string); std::vector<string> P3 = split(P3string, ','); std::vector<size_t> P3pos;
    string P2string; getline(*setsFile, P2string); std::vector<string> P2 = split(P2string, ','); std::vector<size_t> P2pos;
    string P1string; getline(*setsFile, P1string); std::vector<string> P1 = split(P1string, ','); std::vector<size_t> P1pos;
    
    // Now go through the vcf and calculate D
    int AABA = 0; int XXAA = 0; int BBBA = 0; int indels = 0; // int ABBABABA = 0;
    int totalVariantNumber = 0; int usedVariantsCounter = 0;
    double Dnumerator = 0; double Ddenominator = 0;
    double lastVarsDnum = 0; double lastVarsDdenom = 0;
    int lastPrint = 0;
    std::vector<double> regionDs;
    std::vector<string> sampleNames;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            if (!opt::bAaEqO) { Opos = locateSet(sampleNames, outgroup); }
            P3pos = locateSet(sampleNames, P3);
            P2pos = locateSet(sampleNames, P2); P1pos = locateSet(sampleNames, P1);
            
            
            if (!opt::bAaEqO) { std::cerr << "Outgroup: "; print_vector_stream(outgroup, std::cerr); } else { std::cerr << "Outgroup: " << outgroupString << std::endl; }
            std::cerr << "P3: "; print_vector_stream(P3, std::cerr);
            std::cerr << "P2: "; print_vector_stream(P2, std::cerr);
            std::cerr << "P1: "; print_vector_stream(P1, std::cerr);
        } else {
            totalVariantNumber++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (info[0] != "INDEL") {
                string AA = split(info[info.size()-1],'=')[1];
                if (!opt::bAaEqO) {
                    FourSetCounts c;
                    if (AA == fields[3]) {
                        c = getFourSetVariantCounts(fields,P1pos,P2pos,P3pos,Opos,"ref");
                    } else if (AA == fields[4]) {
                        c = getFourSetVariantCounts(fields,P1pos,P2pos,P3pos,Opos,"alt");
                    }
                    Dnumerator += ((1-c.set1daAF)*c.set2daAF*c.set3daAF*(1-c.set4daAF)) - (c.set1daAF*(1-c.set2daAF)*c.set3daAF*(1-c.set4daAF));
                    Ddenominator += ((1-c.set1daAF)*c.set2daAF*c.set3daAF*(1-c.set4daAF)) + (c.set1daAF*(1-c.set2daAF)*c.set3daAF*(1-c.set4daAF));
                } else {
                    ThreeSetCounts c;
                    if (AA == fields[3]) {
                        c = getThreeSetVariantCounts(fields,P1pos,P2pos,P3pos,"ref");
                    } else if (AA == fields[4]) {
                        c = getThreeSetVariantCounts(fields,P1pos,P2pos,P3pos,"alt");
                    }
                    if (c.set3daAF == 0) {
                        XXAA++; continue;
                    } else if (c.set1daAF == 0 && c.set2daAF == 0) {
                        AABA++; continue;
                    } else if (c.set1daAF == 1 && c.set2daAF == 1 && c.set3daAF == 1) {
                        BBBA++; continue;
                    } else {
                        usedVariantsCounter++;
                        // Green et al. (2010) eq. S15.2
                        Dnumerator += ((1-c.set1daAF)*c.set2daAF*c.set3daAF) - (c.set1daAF*(1-c.set2daAF)*c.set3daAF);
                        Ddenominator += ((1-c.set1daAF)*c.set2daAF*c.set3daAF) + (c.set1daAF*(1-c.set2daAF)*c.set3daAF);
                        lastVarsDnum += ((1-c.set1daAF)*c.set2daAF*c.set3daAF) - (c.set1daAF*(1-c.set2daAF)*c.set3daAF);
                        lastVarsDdenom += ((1-c.set1daAF)*c.set2daAF*c.set3daAF) + (c.set1daAF*(1-c.set2daAF)*c.set3daAF);
                        if (Dnumerator == 0) {
                            std::cerr << "P1:" << c.set1daAF << " P2:" << c.set2daAF << " P3:" << c.set3daAF << std::endl;
                        }
                        assert(Dnumerator != 0);
                    }
                }
                // if (totalVariantNumber % 100000 == 0) { std::cerr << Dnumerator << std::endl; }
            } else {
                indels++;
            }
            if (usedVariantsCounter % 5000 == 0 && usedVariantsCounter != lastPrint) {
            //if (totalVariantNumber % 100000 == 0) {
                assert(XXAA + AABA + BBBA + indels + usedVariantsCounter == totalVariantNumber);
                if (usedVariantsCounter > 30000) {
                    double Dstd_err = jackknive_std_err(regionDs);
                    std::cerr << totalVariantNumber << " variants processed. " << usedVariantsCounter << " variants used. \tD=" << (double)Dnumerator/Ddenominator << " std_err=" <<Dstd_err << std::endl;
                } else {
                    std::cerr << totalVariantNumber << " variants processed. " << usedVariantsCounter << " variants used. \tD=" << (double)Dnumerator/Ddenominator << std::endl;
                }
                std::cerr << "Last used 5000 variants \t\t\t\tD=" << lastVarsDnum/lastVarsDdenom << std::endl;
                // std::cerr << "AAAA=" << XXAA << "; AABA=" << AABA << "; BBBA=" << BBBA << std::endl;
                regionDs.push_back(lastVarsDnum/lastVarsDdenom);
                lastVarsDnum = 0; lastVarsDdenom = 0;
                lastPrint = usedVariantsCounter;
            }
        }
    }
    
    double Dstd_err = jackknive_std_err(regionDs);
    std::cerr << std::endl;
    std::cerr << totalVariantNumber << " variants processed. D=" << (double)Dnumerator/Ddenominator << " std_err=" << Dstd_err << std::endl;
    
    
}


int abbaBabaMain(int argc, char** argv) {
    parseAbbaBabaOptions(argc, argv);
    
    if (opt::bFrequency) { frequencyData(); }
    else { singleSequence(); }
    
    return 0;
    
}

void parseAbbaBabaOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 's': arg >> opt::sampleNameFile; break;
            case 'f': opt::bFrequency = true; break;
            case OPT_AA_EQ_O: opt::bAaEqO = true; break;
            case 'h':
                std::cout << ABBA_USAGE_MESSAGE;
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
        std::cout << "\n" << ABBA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}



