//
//  evo_shared_variation.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 25/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "evo_shared_variation.h"
#include "process_vcf_stats_utils.h"

#define SUBPROGRAM "sharedVariation"

static const char *SHAREDVAR_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE.vcf.gz SETS.txt\n"
"Search for shared polymorphic sites between groups, and also for shared heterozyhous sites\n"
"The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID\n"
"\n"
"       -h, --help                                  display this help and exit\n"
"       -l, --sharedVarLocations=S1,S2              (optional) should be a pair of populations for which the location of shared\n"
"                                                   polymorphic sites will be written into a file\n"
"       -n, --run-name                              run-name will be included in the output file name\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "hn:l:";

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "sharedVarLocations",   required_argument, NULL, 'l' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static std::vector<string> getLocsFor;
    static string runName = "";
}

int sharedVarMain(int argc, char** argv) {
    parseSharedVarOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    string fileRoot = stripExtension(opt::setsFile);
    std::cerr << "Looking shared variation between samples in: " << opt::vcfFile << std::endl;
    std::cerr << "The groups are defined in: " << opt::setsFile << std::endl;
    
    std::map<string, std::vector<string>> speciesToIDsMap;
    std::map<string, string> IDsToSpeciesMap;
    std::map<string, std::vector<size_t>> speciesToPosMap;
    std::map<size_t, string> posToSpeciesMap;
    
    // Get the sample sets
    while (getline(*setsFile, line)) {
        std::vector<string> ID_Species = split(line, '\t');
        speciesToIDsMap[ID_Species[1]].push_back(ID_Species[0]);
        IDsToSpeciesMap[ID_Species[0]] = ID_Species[1];
    }
    
    // Get a vector of set names (usually species)
    std::vector<string> species;
    for(std::map<string,std::vector<string>>::iterator it = speciesToIDsMap.begin(); it != speciesToIDsMap.end(); ++it) {
        if ((it->first) != "Outgroup" && it->first != "xxx") {
            species.push_back(it->first);
            // std::cerr << it->first << std::endl;
        }
    } std::cerr << "There are " << species.size() << " sets (excluding the Outgroup)" << std::endl;
    
    std::ofstream* sharedLocsFile;
    if (opt::getLocsFor.size() == 2) {
        sharedLocsFile = new std::ofstream(opt::runName+"sharedVariationLocation_" + opt::getLocsFor[0] + "_" + opt::getLocsFor[1] + ".txt");
    }

    string sharedPerIndividualN = opt::runName + "sharedHets_perIndividual.txt";
    string sharedBetweenGroupsN = "sharedVariationBetween_" + fileRoot + "_" + opt::runName + ".txt";
    string sharedPerIndividualNscaled = opt::runName + "sharedHets_perIndividual_scaled.txt";
    string sharedBetweenGroupsNscaled = "sharedVariationBetween_" + fileRoot + "_" + opt::runName + "_scaled.txt";
    std::ofstream* sharedPerIndividual = new std::ofstream(sharedPerIndividualN.c_str());
    std::ofstream* sharedBetweenGroups = new std::ofstream(sharedBetweenGroupsN.c_str());
    std::ofstream* sharedPerIndividualScaled = new std::ofstream(sharedPerIndividualNscaled.c_str());
    std::ofstream* sharedBetweenGroupsScaled = new std::ofstream(sharedBetweenGroupsNscaled.c_str());
    std::vector<std::vector<double> > hetMatrix; std::vector<std::vector<double> > hetMatrixMissing;
    std::vector<std::vector<double> > sharedBetweenGroupsMatrix; std::vector<std::vector<double> > sharedBetweenGroupsMatrixMissing;
    
    int totalVariantNumber = 0; int reportProgressEvery = 10000;
    std::vector<string> sampleNames; int numSamples = 0;
    std::vector<string> fields;
    std::map<std::string, double> loc_pval;
    std::clock_t start; double durationOverall;
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sN(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            sampleNames = sN;
            for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
                posToSpeciesMap[i] = IDsToSpeciesMap[sampleNames[i]];
            }
            numSamples = (int)sampleNames.size();
            
            // Iterate over all the keys in the map to find the samples in the VCF:
            // Give an error if no sample is found for a species:
            for(std::map<string, std::vector<string>>::iterator it = speciesToIDsMap.begin(); it != speciesToIDsMap.end(); ++it) {
                string sp =  it->first;
                //std::cerr << "sp " << sp << std::endl;
                std::vector<string> IDs = it->second;
                std::vector<size_t> spPos = locateSet(sampleNames, IDs);
                if (spPos.empty()) {
                    std::cerr << "Did not find any samples in the VCF for \"" << sp << "\"" << std::endl;
                    assert(!spPos.empty());
                }
                speciesToPosMap[sp] = spPos;
            }
            
            initialize_matrix_double(hetMatrix, numSamples); initialize_matrix_double(hetMatrixMissing, numSamples);
            initialize_matrix_double(sharedBetweenGroupsMatrix, (int)species.size());
            initialize_matrix_double(sharedBetweenGroupsMatrixMissing, (int)species.size());
            start = std::clock();
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants in " << durationOverall << "secs" << std::endl;
            }
            fields = split(line, '\t');
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());

            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1) {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            GeneralSetCounts* c = new GeneralSetCounts(speciesToPosMap, (int)genotypes.size());
            c->getSetVariantCounts(genotypes, posToSpeciesMap);
            genotypes.clear(); genotypes.shrink_to_fit();
            
            // Put the number of hets in each sample on the diagonal
            for (int i = 0; i < numSamples; i++) {
                if (c->individualsWithVariant[i] == 1)
                    hetMatrix[i][i]++;
                else if (c->individualsWithVariant[i] == -1) // Missing data for this individual
                    hetMatrixMissing[i][i]++;
            }
            // Now look if the het is shared between individuals
            for (int i = 0; i < (numSamples-1); i++) {
                int iAllele = c->individualsWithVariant[i];
                if (iAllele == 1) {
                    for (int j = i+1; j != numSamples; j++) {
                        if (c->individualsWithVariant[j] == 1)
                            hetMatrix[j][i]++;
                        else if (c->individualsWithVariant[j] == -1)
                            hetMatrixMissing[j][i]++;
                    }
                } else if (iAllele == -1) {
                    for (int j = i+1; j != numSamples; j++) {
                        hetMatrixMissing[j][i]++;
                    }
                }
            }
            
            // Now look at the numbers a polymorphic sites in groups and shared between groups:
            std::vector<double> allPs;
            for (int i = 0; i < (int)species.size(); i++) {
                allPs.push_back(c->setAAFs.at(species[i]));
                double p_Si = allPs[i];
                if (p_Si > 0 && p_Si < 1)
                    sharedBetweenGroupsMatrix[i][i]++;
                else if (p_Si == -1)
                    sharedBetweenGroupsMatrixMissing[i][i]++;
            }
            for (int i = 0; i < ((int)species.size()-1); i++) {
                double p_Si = allPs[i];
                if (p_Si > 0 && p_Si < 1) {
                    for (int j = i+1; j != (int)species.size(); j++) {
                        double p_Sj = allPs[j];
                        if (p_Sj > 0 && p_Sj < 1) {
                            sharedBetweenGroupsMatrix[j][i]++;
                            if (opt::getLocsFor.size() == 2) {
                                if ((species[i] == opt::getLocsFor[0] && species[j] == opt::getLocsFor[1]) || (species[i] == opt::getLocsFor[1] && species[j] == opt::getLocsFor[0]))
                                    *sharedLocsFile << fields[0] << "\t" << fields[1];
                            }
                        }
                        else if (p_Sj == -1)
                            sharedBetweenGroupsMatrixMissing[j][i]++;
                    }
                } else if (p_Si == -1) { // Entirely missing data for this species/population
                    for (int j = i+1; j != (int)species.size(); j++) {
                        sharedBetweenGroupsMatrixMissing[j][i]++;
                    }
                }
            }
            delete c;
        }
    }
    
    // print het sharing per individual
    print_vector(sampleNames, *sharedPerIndividual);
    print_matrix<const std::vector<std::vector<double> >&>(hetMatrix, *sharedPerIndividual);
    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numSamples; j++) {
            double proportionUsed = 1 - (hetMatrixMissing[j][i]/totalVariantNumber);
            hetMatrix[j][i] = hetMatrix[j][i]/proportionUsed;
        }
    }
    print_vector(sampleNames, *sharedPerIndividualScaled);
    print_matrix<const std::vector<std::vector<double> >&>(hetMatrix, *sharedPerIndividualScaled);
    
    // Print variant sharing between groups
    print_vector(species, *sharedBetweenGroups);
    print_matrix<const std::vector<std::vector<double> >&>(sharedBetweenGroupsMatrix, *sharedBetweenGroups);
    for (int i = 0; i < (int)species.size(); i++) {
        for (int j = 0; j < (int)species.size(); j++) {
            double proportionUsed = 1 - (sharedBetweenGroupsMatrixMissing[j][i]/totalVariantNumber);
            sharedBetweenGroupsMatrix[j][i] = sharedBetweenGroupsMatrix[j][i]/proportionUsed;
        }
    }
    print_vector(species, *sharedBetweenGroupsScaled);
    print_matrix<const std::vector<std::vector<double> >&>(sharedBetweenGroupsMatrix, *sharedBetweenGroupsScaled);
    
    
    return 0;
}

void parseSharedVarOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'n': arg >> opt::runName; break;
            case 'l': opt::getLocsFor = split(arg.str(), ','); break;
            case '?': die = true; break;
            case 'h': std::cout << SHAREDVAR_USAGE_MESSAGE;
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
        std::cout << "\n" << SHAREDVAR_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}
