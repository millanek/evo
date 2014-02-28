//
//  process_vcf_search_sex.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 25/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_search_sex.h"
#include "process_vcf_stats_utils.h"

#define SUBPROGRAM "sex-search"

static const char *SEX_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] FILTERED_VCF_FILE SAMPLE_SETS.txt\n"
"Looking for regions of the genome corresponding to the X and Y sex determining sequences\n"
"based on differences in coverage as recorded in over variants in the multi-sample vcf file\n"
"The SAMPLE_SETS.txt file should have exactly two lines; first line with all the males and second line with all the females\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -w, --window NUM                        width of the window in NUM variants (default 200)\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the VCF file\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "hs:w:";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "window",   required_argument, NULL, 'w' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genderFile;
    static string sampleNameFile;
    static int window = 200;
}

SetDepths getSetVariantDepths(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci) {
    SetDepths thisVariantDepths;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        std::vector<string> genotype_field = split(fields[i], ':');
        if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantDepths.set1Depths.push_back(atoi(genotype_field[2].c_str())); }
        if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantDepths.set2Depths.push_back(atoi(genotype_field[2].c_str())); }
    }
    thisVariantDepths.fillMeanDepths();
    return thisVariantDepths;
}

 

int sexSearchMain(int argc, char** argv) {
    parseSexSearchOptions(argc, argv);
    string fileRoot = stripExtension(opt::vcfFile);
    
    std::cerr << "Looking for the X and Y sex determining sequences using data from: " << opt::vcfFile << std::endl;
    std::cerr << "and male/female status defined in: " << opt::genderFile << std::endl;
    std::cerr << "The window width is: " << opt::window << " variants" << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* genderFile = new std::ifstream(opt::genderFile.c_str());
    
    string pValFileName = "sex_t_pvals_w" + numToString(opt::window) + ".txt";
    std::ofstream* pValFile = new std::ofstream(pValFileName.c_str());
    
    string maleString; string femaleString;
    getline(*genderFile, maleString);
    getline(*genderFile, femaleString);
    std::vector<string> males = split(maleString, ',');
    std::vector<string> females = split(femaleString, ',');
    std::sort(males.begin(),males.end());
    std::sort(females.begin(),females.end());
    
    
    // Prepare an object for output files
    
    string line;
    string currentScaffoldNum = "";
    string currentScaffoldReference;
    std::vector<string> sampleNames;
    string thisScaffoldName;
    std::vector<size_t> maleLoci;
    std::vector<size_t> femaleLoci;
    std::vector<SetDepths> setDepthsVector;
    std::vector<double> set1MeanDepths;
    std::vector<double> set2MeanDepths;
    unsigned int processedVariantCounter = 0;
    unsigned int regionStart = 1;
    std::vector<int> YlikeCounts;
    std::vector<int> nonZeroYlikeCounts;
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            continue;
        } else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            maleLoci = locateSet(sampleNames, males);
            femaleLoci = locateSet(sampleNames, females);
        } else {
            processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            
            SetDepths depths = getSetVariantDepths(fields,maleLoci,femaleLoci);
            setDepthsVector.push_back(depths);
            set1MeanDepths.push_back(depths.set1_mean_depth);
            set2MeanDepths.push_back(depths.set2_mean_depth);
            
            if (processedVariantCounter % opt::window == 0) {
                int Ylike = 0;
                for (std::vector<SetDepths>::size_type i = 0; i != setDepthsVector.size(); i++) {
                    // Could this be a Y region
                    if(setDepthsVector[i].set1_mean_depth > 2 && setDepthsVector[i].set2_mean_depth <= 1) {
                        Ylike++;
                    }
                }
                
                // std::cerr << "Region: " << fields[0] << "\t" << regionStart << "-" << fields[1] << "processed" << std::endl;
                // std::cerr << Ylike << " Ylike variants" << std::endl;
                YlikeCounts.push_back(Ylike);
                if (Ylike > (opt::window/2)) {
                    std::cout << "Putative Y region: " << fields[0] << "\t" << regionStart << "-" << fields[1] << std::endl;
                }
                if (Ylike > 0) {
                    nonZeroYlikeCounts.push_back(Ylike);
                    std::cerr << "Region: " << fields[0] << "\t" << regionStart << "-" << fields[1] << "processed" << std::endl;
                    std::cerr << Ylike << " Ylike variants" << std::endl;
                }
                
                
               /* if (regionStart == 182909) {
                   // std::cerr << "Here " << set1MeanDepths.size() << "\t" << set1MeanDepths[50] << std::endl;
                  //  print_vector_stream(set1MeanDepths, std::cerr, '\n');
                  //  std::cerr << "Here " << set2MeanDepths.size() << std::endl;
                  //  print_vector_stream(depths.set1Depths, std::cerr, '\n');
                    //std::cerr << "Here " << set2MeanDepths.size() << std::endl;
                    //print_vector_stream(set2MeanDepths, std::cerr, '\n');
                    //std::cerr << "Here " << set2MeanDepths.size() << std::endl;
                    //print_vector_stream(depths.set2Depths, std::cerr, '\n');
                }
                */
                // alternative hypothesis: true difference in means is greater than 5
                double pval_t = two_sample_t(set1MeanDepths,set2MeanDepths);
                
                if (pval_t < 0.001) {
                    string loc = fields[0] + "\t" + numToString(regionStart) + "\t" + fields[1];
                    *pValFile << loc << "\t" << pval_t << "\t" << vector_average(set1MeanDepths) << "\t" << vector_average(set2MeanDepths) << std::endl;
                }
                
                //std::cerr << "Region: " << fields[0] << "\t" << regionStart << "-" << fields[1] << "processed" << std::endl;
                regionStart = atoi(fields[1].c_str());
                setDepthsVector.clear();
                set1MeanDepths.clear(); set2MeanDepths.clear();
                
            }
            if (processedVariantCounter % 100000 == 0)
                std::cerr << processedVariantCounter << " variants processed" << std::endl;
        }
    }
    
    std::ofstream* YlikeCountsDistribution = new std::ofstream("YlikeCounts.txt");
    print_vector(YlikeCounts, *YlikeCountsDistribution, '\n');
    std::ofstream* nonZeroYlikeCountsDistribution = new std::ofstream("nonZeroYlikeCounts.txt");
    print_vector(YlikeCounts, *nonZeroYlikeCountsDistribution, '\n');
    return 0;
}

void parseSexSearchOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 's': arg >> opt::sampleNameFile; break;
            case 'w': arg >> opt::window; break;
            case '?': die = true; break;
            case 'h': std::cout << SEX_USAGE_MESSAGE;
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
        std::cout << "\n" << SEX_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::genderFile = argv[optind++];
}
