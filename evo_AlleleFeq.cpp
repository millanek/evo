//
//  evo_AlleleFeq.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 26/05/2020.
//  Copyright © 2020 Milan Malinsky. All rights reserved.
//

#include "evo_AlleleFeq.h"
#include "process_vcf_annotation_tools.h"
#include <deque>

#define SUBPROGRAM "alleleFreq"

#define DEBUG 1

static const char *AF_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt\n"
"Calculate the Allele Frequencies per population/species from a VCF \n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -g, --use-genotype-probabilities        (optional) use probabilities (GP tag) or calculate them from likelihoods (GL or PL tags) using a Hardy-Weinberg prior\n"
"                                               the probabilities are used to estimate allele frequencies in each population/species\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:";


static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "use-genotype-probabilities", no_argument, NULL, 'g'},
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string runName = "out";
    static bool useGenotypeProbabilities = false;
}


int AFmain(int argc, char** argv) {
    parseAFoptions(argc, argv);
    string line; // for reading the input files
    

    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    
    std::map<string, std::vector<string>> popToIDsMap;
    std::map<string, string> IDsToPopMap;
    std::map<string, std::vector<size_t>> popToPosMap;
    std::map<size_t, string> posToPopMap;
    
    // Get the sample sets
    while (getline(*setsFile, line)) {
        // std::cerr << line << std::endl;
        std::vector<string> ID_Pop = split(line, '\t');
        popToIDsMap[ID_Pop[1]].push_back(ID_Pop[0]);
        IDsToPopMap[ID_Pop[0]] = ID_Pop[1];
        //std::cerr << ID_Species[1] << "\t" << ID_Species[0] << std::endl;
    }
    // Get a vector of set names (usually populations/species)
    std::vector<string> populations;
    for(std::map<string,std::vector<string>>::iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
        populations.push_back(it->first);
        // std::cerr << it->first << std::endl;
    } std::cerr << "There are " << populations.size() << " populations " << std::endl;

    int totalVariantNumber = 0;
    std::vector<string> sampleNames; std::vector<std::string> fields;
    int reportProgressEvery = 10000; string chr; string coord; double coordDouble;
    std::clock_t start; std::clock_t startGettingCounts;
    double durationOverall; double durationGettingCounts; 
    //std::ofstream* outFileAF = new std::ofstream(stripExtension(opt::setsFile) + "_" + opt::runName + "_AF" + ".txt");
    std::ostream* outFileAF = createWriter(stripExtension(opt::setsFile) + "_" + opt::runName + "_AF" + ".txt.gz");
    
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            // print_vector_stream(sampleNames, std::cerr);
            for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
                posToPopMap[i] = IDsToPopMap[sampleNames[i]];
            }
            // Iterate over all the keys in the map to find the samples in the VCF:
            // Give an error if no sample is found for a species:
            for(std::map<string, std::vector<string>>::iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
                string sp =  it->first;
                //std::cerr << "sp " << sp << std::endl;
                std::vector<string> IDs = it->second;
                std::vector<size_t> spPos = locateSet(sampleNames, IDs);
                if (spPos.empty()) {
                    std::cerr << "Did not find any samples in the VCF for \"" << sp << "\"" << std::endl;
                    assert(!spPos.empty());
                }
                popToPosMap[sp] = spPos;
            }
            start = std::clock();
            //  std::cerr << " " << std::endl;
            //  std::cerr << "Outgroup at pos: "; print_vector_stream(speciesToPosMap["Outgroup"], std::cerr);
            //  std::cerr << "telvit at pos: "; print_vector_stream(speciesToPosMap["telvit"], std::cerr);
        } else {
            totalVariantNumber++;
            if (totalVariantNumber == 1) {
                *outFileAF << "chr" << "\t" << "coord" << "\t" << "ref" << "\t" << "alt";
                for(std::map<string,std::vector<size_t>>::iterator iter =  popToPosMap.begin(); iter != popToPosMap.end(); ++iter) {
                    *outFileAF << "\t" << iter->first;
                }
                *outFileAF << "\n";
            }
            if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants in " << durationOverall << "secs" << std::endl;
            }
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1]; coordDouble = stringToDouble(coord);
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1 || altAllele == "*") {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            startGettingCounts = std::clock();
            GeneralSetCounts* c = new GeneralSetCounts(popToPosMap, (int)genotypes.size());
            c->getSetVariantCountsSimple(genotypes, posToPopMap);
            durationGettingCounts = ( std::clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            // std::cerr << "Here:" << totalVariantNumber << std::endl;

            
            *outFileAF << chr << "\t" << coord << "\t" << refAllele << "\t" << altAllele;
            if (opt::useGenotypeProbabilities) {
                int likelihoodsOrProbabilitiesTagPosition = c->checkForGenotypeLikelihoodsOrProbabilities(fields);
                if (likelihoodsOrProbabilitiesTagPosition == LikelihoodsProbabilitiesAbsent) {
                    printMissingLikelihoodsWarning(fields[0], fields[1]);
                    opt::useGenotypeProbabilities = false;
                } else c->getAFsFromGenotypeLikelihoodsOrProbabilities(genotypes,posToPopMap,likelihoodsOrProbabilitiesTagPosition);
                for(std::map<string,double>::iterator iter =  c->setAAFsFromLikelihoods.begin(); iter != c->setAAFsFromLikelihoods.end(); ++iter) {
                    *outFileAF << "\t" << iter->second;
                }
            
            } else {
                for(std::map<string,double>::iterator iter =  c->setAAFs.begin(); iter != c->setAAFs.end(); ++iter) {
                    *outFileAF << "\t" << iter->second;
                }
            }
            *outFileAF << "\n";

            genotypes.clear(); genotypes.shrink_to_fit();
            
            delete c;
        }
    }
    
    outFileAF->flush();
    return 0;
}



void parseAFoptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'g': opt::useGenotypeProbabilities = true; break;
            case 'h':
                std::cout << AF_USAGE_MESSAGE;
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
        std::cout << "\n" << AF_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}

