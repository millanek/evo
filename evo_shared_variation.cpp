//
//  evo_shared_variation.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 25/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "evo_shared_variation.h"
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

int sharedVarMain(int argc, char** argv) {
    parseSharedVarOptions(argc, argv);
    string fileRoot = stripExtension(opt::sampleSets);
    GenericSampleSets* sampleSets;
    std::cerr << "Looking shared variation between samples in: " << opt::vcfFile << std::endl;
    std::cerr << "The groups are defined in: " << opt::sampleSets << std::endl;
    
    // Open connection to read from the vcf file
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::sampleSets.c_str());
    string sharedPerIndividualN = opt::runName + "sharedHets_perIndividual.txt";
    string sharedBetweenGroupsN = "sharedVariationBetween_" + fileRoot + "_" + opt::runName + ".txt";
    string riverineSharedWithLakesN = "riverineSharedWithLakes_" + fileRoot + "_" + opt::runName + ".txt";
    std::ofstream* sharedPerIndividual = new std::ofstream(sharedPerIndividualN.c_str());
    std::ofstream* sharedBetweenGroups = new std::ofstream(sharedBetweenGroupsN.c_str());
    std::ofstream* riverineSharedWithLakes;
    std::vector<std::vector<double> > hetMatrix;
    std::vector<std::vector<double> > sharedBetweenGroupsMatrix;
    std::vector<std::vector<double> > riverineWithLakesMatrix;
    
    int numChromosomes;
    int totalVariantNumber = 0;
    string line;
    std::vector<string> sampleNames; int numSamples = 0;
    std::vector<string> fields;
    std::map<std::string, double> loc_pval;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
 
        } else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            numSamples = (int)fields.size() - NUM_NON_GENOTYPE_COLUMNS;
            numChromosomes = numSamples * 2;
            //std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
            
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            
            // std::cerr << "here: " << std::endl;
            sampleSets = new GenericSampleSets(setsFile,sampleNames);
            for (int i = 0; i < (int)sampleSets->sets.size(); i++) {
                std::cerr << "Set " << sampleSets->allSetNames[i] << " loci: " << std::endl;
                print_vector_stream(sampleSets->sets[sampleSets->allSetNames[i]], std::cerr);
            }
            //std::cerr << "Set " << (++sampleSets->sets.begin())->first << " loci: " << std::endl;
            //print_vector_stream((++sampleSets->sets.begin())->second, std::cerr);
            initialize_matrix_double(hetMatrix, numSamples);
            initialize_matrix_double(sharedBetweenGroupsMatrix, (int)sampleSets->sets.size());
            if (sampleSets->sets.count("Riverine") == 1) {
                riverineSharedWithLakes = new std::ofstream(riverineSharedWithLakesN.c_str());
                initialize_matrix_double(riverineWithLakesMatrix, (int)sampleSets->sets.size()-1, (int)sampleSets->sets["Riverine"].size());
            }
            
        } else {
            totalVariantNumber++;
            
            std::vector<std::string> fields = split(line, '\t');
            ManySetCountsGeneric setCounts = getGenericSetVariantCounts(fields, *sampleSets);
            
            // Put the number of hets in each sample on the diagonal
            for (int i = 0; i < numSamples; i++) {
                if (setCounts.individualsWithVariant[i] == 1)
                    hetMatrix[i][i]++;
            }
            // Now look if the het is shared between individuals
            for (int i = 0; i < (numSamples-1); i++) {
                for (int j = i+1; j != numSamples; j++) {
                    if (setCounts.individualsWithVariant[i] == 1 && setCounts.individualsWithVariant[j] == 1)
                        hetMatrix[j][i]++;
                }
            }
            
            // Now look at the numbers a polymorphic sites in groups and shared between groups:
            for (int i = 0; i < (int)sampleSets->sets.size(); i++) {
                if (setCounts.allSetCounts[sampleSets->allSetNames[i]].isPolymorhic())
                    sharedBetweenGroupsMatrix[i][i]++;
            }
            for (int i = 0; i < ((int)sampleSets->sets.size()-1); i++) {
                for (int j = i+1; j != (int)sampleSets->sets.size(); j++) {
                    if (setCounts.allSetCounts[sampleSets->allSetNames[i]].isPolymorhic() && setCounts.allSetCounts[sampleSets->allSetNames[j]].isPolymorhic())
                    sharedBetweenGroupsMatrix[j][i]++;
                }
            }
            //std::cerr << "Hello: ";
            
            if (sampleSets->sets.count("Riverine") == 1) {
                // And finally, sharing between riverines and the great lakes:
                for (int i = 0; i < (int)sampleSets->sets["Riverine"].size(); i++) {
                    //std::cerr << "Hello: ";
                    // If this riverine individual is a het:
                    if (setCounts.individualsWithVariant[sampleSets->sets["Riverine"][i]] == 1) {
                        //std::cerr << "\tHet\t";
                        riverineWithLakesMatrix[i][3]++;
                        if (setCounts.allSetCounts["Malawi"].isPolymorhic())
                            riverineWithLakesMatrix[i][0]++;
                        if (setCounts.allSetCounts["Tanganyika"].isPolymorhic())
                            riverineWithLakesMatrix[i][1]++;
                        if (setCounts.allSetCounts["Victoria"].isPolymorhic())
                            riverineWithLakesMatrix[i][2]++;
                    }
                   // std::cerr << "\tBye: " << std::endl;
                }
            }
            
            
            if (totalVariantNumber % 100000 == 0)
                std::cerr << "Processed " << totalVariantNumber << " variants" << std::endl;
        }
    }
    
    // print het sharing per individual
    print_vector(sampleNames, *sharedPerIndividual);
    print_matrix<const std::vector<std::vector<double> >&>(hetMatrix, *sharedPerIndividual);
    
    // Print variant sharing between groups
    print_vector(sampleSets->allSetNames, *sharedBetweenGroups);
    print_matrix<const std::vector<std::vector<double> >&>(sharedBetweenGroupsMatrix, *sharedBetweenGroups);
    
    // Print variant sharing between riverine individuals and each of the great lakes:
    if (sampleSets->sets.count("Riverine") == 1) {
        *riverineSharedWithLakes << "Malawi\tTanganyika\tVictoria\tTotalHets" << std::endl;
        for (int i = 0; i < (int)sampleSets->sets["Riverine"].size(); i++) {
            *riverineSharedWithLakes << sampleNames[sampleSets->sets["Riverine"][i]] << "\t"; print_vector(riverineWithLakesMatrix[i], *riverineSharedWithLakes);
        }
    }
    return 0;
}

void parseSharedVarOptions(int argc, char** argv) {
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
