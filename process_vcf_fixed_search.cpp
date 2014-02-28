//
//  process_vcf_fixed_search.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 25/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_fixed_search.h"
#include "process_vcf_stats_utils.h"

#define SUBPROGRAM "fixed-search"

static const char *FIXED_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] FILTERED_VCF_FILE SAMPLE_SETS.txt\n"
"Search for fixed (or nearly fixed) differences between two sets of samples\n"
"This will be much faster if singletons (and probably doubletons) are filtered out from the VCF file before running this\n"
"The SAMPLE_SETS.txt file should have exactly two lines, with each line defining one of the two sample sets\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -p, --pval-output-cutoff                maximum p-value for output in the _pvals.txt file; (default 0.01)\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the VCF file\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "hn:s:p:";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "pval-output-cutoff",   required_argument, NULL, 'p' },
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
    static double pvalCutoff = 0.01;
}

SetCounts getSetVariantCounts(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci) {
    SetCounts thisVariantCounts;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1Count++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2Count++; }
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1Count++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2Count++; }
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
    }
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    int set1WithoutVariant = (int)(set1_loci.size() * 2) - thisVariantCounts.set1Count;
    int set2WithoutVariant = (int)(set2_loci.size() * 2) - thisVariantCounts.set2Count;
    /*
    if (fields[1] == "2520") {
        std::cerr << thisVariantCounts.set1Count << "\t" << set1WithoutVariant << "\t" << thisVariantCounts.set2Count << "\t" << set2WithoutVariant << std::endl;
    }
     */
    if ((thisVariantCounts.set1Count != 0 || thisVariantCounts.set2Count != 0) && (set1WithoutVariant != 0 || set2WithoutVariant != 0)) {
        if ((set1_loci.size() * 2) + (set2_loci.size() * 2) <= 60) {
            thisVariantCounts.fisher_pval = fisher_exact(thisVariantCounts.set1Count,set1WithoutVariant , thisVariantCounts.set2Count, set2WithoutVariant);
            thisVariantCounts.chi_sq_pval = pearson_chi_sq_indep(thisVariantCounts.set1Count,set1WithoutVariant , thisVariantCounts.set2Count, set2WithoutVariant);
        } else {
            thisVariantCounts.chi_sq_pval = pearson_chi_sq_indep(thisVariantCounts.set1Count,set1WithoutVariant , thisVariantCounts.set2Count, set2WithoutVariant);
        }
    }
    return thisVariantCounts;
}

int getNumHets(SetCounts& counts) {
    int num_hets = 0;
    for (std::vector<std::vector<int> >::size_type i = 0; i < counts.individualsWithVariant.size(); i++) {
        if (counts.individualsWithVariant[i] == 1)
            num_hets++;
    }
    return num_hets;
}

int fixedSearchMain(int argc, char** argv) {
    parseFixedSearchOptions(argc, argv);
    string fileRoot = stripExtension(opt::sampleSets);
    
    std::cerr << "Looking for fixed variants in: " << opt::vcfFile << std::endl;
    std::cerr << "Between two sets of samples defined in: " << opt::sampleSets << std::endl;
    
    // Open connection to read from the vcf file
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::sampleSets.c_str());
    string VarFixedFileName = fileRoot + "_" + opt::runName + "_fixed.vcf";
    string VarTier2FileName = fileRoot + "_" + opt::runName + "_tier2.vcf";
    string VarTier3FileName = fileRoot + "_" + opt::runName + "_tier3.vcf";
    string VarAllTiersFileName = fileRoot + "_" + opt::runName + "_all_tiers.vcf";
    string PvalFileName = fileRoot + "_" + opt::runName + "_pvals.txt";
    string misclasFileName = fileRoot + "_" + opt::runName + "_misfits.txt";
    std::ofstream* pFixed = new std::ofstream(VarFixedFileName.c_str());
    std::ofstream* pTier2 = new std::ofstream(VarTier2FileName.c_str());
    std::ofstream* pTier3 = new std::ofstream(VarTier3FileName.c_str());
    std::ofstream* pAllTiers = new std::ofstream(VarAllTiersFileName.c_str());
    std::ofstream* pValFile = new std::ofstream(PvalFileName.c_str());
    std::ofstream* pMisFile = new std::ofstream(misclasFileName.c_str());
    
    string set1String; string set2String;
    getline(*setsFile, set1String);
    getline(*setsFile, set2String);
    std::vector<string> set1 = split(set1String, ',');
    std::vector<string> set2 = split(set2String, ',');
    std::sort(set1.begin(),set1.end());
    std::sort(set2.begin(),set2.end());
    
    
    int numChromosomes;
    int totalVariantNumber = 0;
    string line;
    std::vector<string> sampleNames;
    std::vector<string> fields;
    std::vector<size_t> set1Loci;
    std::vector<size_t> set2Loci;
    std::map<std::string, double> loc_pval;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            *pFixed << line << std::endl;
            *pTier2 << line << std::endl;
            *pTier3 << line << std::endl;
            *pAllTiers << line << std::endl;
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
            set1Loci = locateSet(sampleNames, set1);
            set2Loci = locateSet(sampleNames, set2);
            
            //print_vector_stream(sampleNames, std::cerr);
            std::cerr << sampleNames[134] << std::endl;
            
            std::cerr << "Set1 loci: " << std::endl;
            print_vector_stream(set1Loci, std::cerr);
            std::cerr << "Set2 loci: " << std::endl;
            print_vector_stream(set2Loci, std::cerr);
            
            *pFixed << line << std::endl;
            *pTier2 << line << std::endl;
            *pTier3 << line << std::endl;
            *pAllTiers << line << std::endl;
        } else {
            totalVariantNumber++;
            
            std::vector<std::string> fields = split(line, '\t');
            
            SetCounts counts = getSetVariantCounts(fields,set1Loci,set2Loci);
            
            if (counts.fisher_pval < opt::pvalCutoff || counts.chi_sq_pval < opt::pvalCutoff) {
                string loc = fields[0] + "\t" + fields[1];
                *pValFile << loc << "\t" << counts.fisher_pval << "\t" << counts.chi_sq_pval << "\t" << counts.set1Count << "\t" << (set1Loci.size()*2)-counts.set1Count << "\t" << counts.set2Count << "\t" << (set2Loci.size()*2)-counts.set2Count << std::endl;
                if (counts.fisher_pval < 0.00000000001) {
                    if (counts.set1Count == 18) {
                        string misfits = "";
                        for (std::vector<int>::size_type i = 0; i != counts.individualsWithVariant.size(); i++) {
                            if (std::find(set2Loci.begin(), set2Loci.end(), i) != set2Loci.end()) {
                                if (counts.individualsWithVariant[i] == 1) {
                                    if (misfits == "")
                                        misfits += sampleNames[i]+"(het)";
                                    else
                                        misfits += "," + sampleNames[i]+"(het)";
                                    
                                }
                                if (counts.individualsWithVariant[i] == 2) {
                                    if (misfits == "")
                                        misfits += sampleNames[i]+"(hom)";
                                    else
                                        misfits += "," + sampleNames[i]+"(hom)";
                                    
                                }
                            }
                        }
                        *pMisFile << loc << "\t" << counts.fisher_pval << "\t" << counts.set1Count << "\t" << (set1Loci.size()*2)-counts.set1Count << "\t" << counts.set2Count << "\t" << (set2Loci.size()*2)-counts.set2Count << "\t" << misfits << std::endl;
                    }
                    if (counts.set1Count == 0) {
                        string misfits = "";
                        for (std::vector<int>::size_type i = 0; i != counts.individualsWithVariant.size(); i++) {
                            if (std::find(set2Loci.begin(), set2Loci.end(), i) != set2Loci.end()) {
                                if (counts.individualsWithVariant[i] == 1) {
                                    if (misfits == "")
                                        misfits += sampleNames[i]+"(het)";
                                    else
                                        misfits += "," + sampleNames[i]+"(het)";
                                    
                                }
                                if (counts.individualsWithVariant[i] == 0) {
                                    if (misfits == "")
                                        misfits += sampleNames[i]+"(hom)";
                                    else
                                        misfits += "," + sampleNames[i]+"(hom)";
                                    
                                }
                            }
                        }
                        *pMisFile << loc << "\t" << counts.fisher_pval << "\t" << counts.set1Count << "\t" << (set1Loci.size()*2)-counts.set1Count << "\t" << counts.set2Count << "\t" << (set2Loci.size()*2)-counts.set2Count << "\t" << misfits << std::endl;
                    }
                }
                
            }
            //loc_pval[loc] = counts.pval;

            if ((counts.set1Count == (set1.size()*2) && counts.set2Count == 0) || (counts.set1Count == 0 && counts.set2Count == (set2.size()*2))) {
                *pFixed << line << std::endl;
                *pAllTiers << line << std::endl;
            }
            if ((counts.set1Count == (set1.size()*2) - 1 && counts.set2Count == 0) || (counts.set1Count == (set1.size()*2) && counts.set2Count == 1) || (counts.set1Count == 0 && counts.set2Count == (set2.size()*2) - 1) || (counts.set1Count == 1 && counts.set2Count == (set2.size()*2))) {
                *pTier2 << line << std::endl;
                *pAllTiers << line << std::endl;
            }
            if ((counts.set1Count == (set1.size()*2) - 2 && counts.set2Count == 0) || (counts.set1Count == (set1.size()*2) && counts.set2Count == 2) ||
              (counts.set1Count == (set1.size()*2) - 1 && counts.set2Count == 1) || (counts.set1Count == (set1.size()*2)-1 && counts.set2Count == 1) ||
                (counts.set1Count == 0 && counts.set2Count == (set2.size()*2)-2) || (counts.set1Count == 2 && counts.set2Count == (set2.size()*2)) ||
                (counts.set1Count == 1 && counts.set2Count == (set2.size()*2)-1) || (counts.set1Count == 1 && counts.set2Count == (set2.size()*2)-1)) {
                *pTier3 << line << std::endl;
                *pAllTiers << line << std::endl;
            }
        }
    }
    
    
    
    return 0;
}

void parseFixedSearchOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 's': arg >> opt::sampleNameFile; break;
            case 'p': arg >> opt::pvalCutoff; break;
            case 'n': arg >> opt::runName; break;
            case '?': die = true; break;
            case 'h': std::cout << FIXED_USAGE_MESSAGE;
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
        std::cout << "\n" << FIXED_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::sampleSets = argv[optind++];
}
