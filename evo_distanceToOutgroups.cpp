//
//  evo_distanceToOutgroups.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 14/02/2022.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#include "evo_distanceToOutgroups.h"

#define SUBPROGRAM "DistOutgroups"

#define DEBUG 1

static const char *DISTOUT_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt OUTGROUPS.txt INGROUPS.txt\n"
"Calculate the genetic distance to given outgroups for populations specifified in:\n"
"The POPULATIONS.txt file; This file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The OUTGROUPS.txt should contain names of the populations which are considered outgroups - one per line\n"
"The INGROUPS.txt should contain names of the populations which are considered ingroups - one per line\n"
"for the outgroups, mean distances to the other populations will be calculated, possibly in windows as specified:\n"
"There can be multiple lines and then the program generates multiple ouput files, named like OUGROUP_DIST_SIZE_STEP.txt\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -f, --fixedW sizeKb                     fixed window size (default: 10kb)\n"
"       -i, --allow-indels-and-multiallelics   (optional) also calculate the PBS score for indels, and multiallelics\n"
"                                               for multiallelics, the first alternate allele is considered\n"
"                                               sites where the ALT allele is simply '*' are ignored\n"
"       -r , --region=start,length              (optional) only process a subset of the VCF file\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       --accessibleGenomeBED=BEDfile.bed           (optional) a bed file specifying the regions of the genome where we could call SNPs\n"
"                                                   the program will calculate the number of accessible bases from this\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:f:i";

enum { OPT_ACC_GEN_BED, OPT_AF  };

static const struct option longopts[] = {
    { "fixedW",   required_argument, NULL, 'f' },
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "allow-indels-and-multiallelics",   no_argument, NULL, 'i' },
    { "accessibleGenomeBED", required_argument, NULL, OPT_ACC_GEN_BED },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string outgroupFile;
    static string ingroupFile;
    static string runName = "";
    static int fixedWindowSize = 10000;
    static bool allowIndels = false;
    static string accesibleGenBedFile;
}

double DxyPerSNPfromAFs(double AF1, double AF2) {
    double dxy = AF1*(1-AF2) + AF2*(1-AF1);
    return dxy;
}


int DistOutMain(int argc, char** argv) {
    parseDistOutOptions(argc, argv);
    string line; // for reading the input files
    
    // 1) Dividing up and assigning the populations
    std::vector<string> populations; // Will hold the names of all populations
    std::vector<string> outgroups; // Will hold the names of the outgroup populations
    std::vector<string> ingroups; // Will hold the names of the ingroup populations
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    std::ifstream* outgroupFile = new std::ifstream(opt::outgroupFile.c_str());
    std::ifstream* ingroupFile = new std::ifstream(opt::ingroupFile.c_str());
    
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
    for(std::map<string,std::vector<string>>::iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
        populations.push_back(it->first);
        // std::cerr << it->first << std::endl;
    } std::cerr << "There are " << populations.size() << " populations " << std::endl;
    
    
    // find the ingroups
    while (getline(*ingroupFile,line)) {
        ingroups.push_back(line);
    }
    
    // Get the outgroups
    std::vector<std::ofstream*> outFilesFixedWindow;
    while (getline(*outgroupFile,line)) {
        // std::cerr << line << std::endl;
        std::ofstream* outFileFixedWindow = new std::ofstream(line + "_DIST_" + opt::runName + "_FW" + numToString(opt::fixedWindowSize) + ".txt");
        *outFileFixedWindow << "chr\twStart\twEnd\tSNPs_used\tSNPs_missing\tAccessibleSizeBP" << "\t"; print_vector(ingroups, *outFileFixedWindow);
        //outFile->setf(std::ios_base::fixed); // Avoid scientific notation in the coordinates
        outFilesFixedWindow.push_back(outFileFixedWindow);
        outgroups.push_back(line);
    }
    
    
    // 2) Loading the accessible genome bad file:
    AccessibleGenome* ag;
    if (!opt::accesibleGenBedFile.empty()) {
       std::ifstream* accessibleGenomeBed = new std::ifstream(opt::accesibleGenBedFile);
       std::cerr << "Loading the accessible genome annotation" << std::endl;
       ag = new AccessibleGenome(accessibleGenomeBed);
       std::cerr << "Done" << std::endl;
   }
     // int numAccessibleBP = ag->getAccessibleBPinRegion(sc, i, i+opt::accessibleGenBedWindow);
    
    
    // 3) Praparing the containers to hold the results
    // 3a) for per-SNP results
    std::vector<std::vector<double>> initVectorFixed(ingroups.size());
    std::vector<std::vector<std::vector<double>>> DxyFixedWindowPerSNP(outgroups.size(),initVectorFixed);
    // 3b) for calculated Dxy values
    std::vector<int> ingroupsMissingInit(ingroups.size(),0);
    std::vector<std::vector<int>> missingDist(outgroups.size(),ingroupsMissingInit);
    
    // 3b) for calculated Dxy values
    std::vector<double> ingroupsResultInit(ingroups.size(),0.0);
    std::vector<std::vector<double>> DxyFixedWindowAveraged(outgroups.size(),ingroupsResultInit);
    
    
    // 4) Going through the VCF
    // 4a) Preparing useful variables
    int currentWindowStart = 0; int currentWindowEnd = currentWindowStart + opt::fixedWindowSize;
    int totalVariantNumber = 0;
    std::vector<int> usedVars(outgroups.size(),0); // Will count the number of used variants for each outgroup
    std::vector<int> missingVars(outgroups.size(),0); // Will count the number of missing variants for each outgroup
    std::vector<string> sampleNames; std::vector<std::string> fields;
    int reportProgressEvery = 1000; string chr; string coord; double coordDouble;
    std::clock_t start; std::clock_t startGettingCounts; std::clock_t startCalculation;
    double durationOverall; double durationGettingCounts; double durationCalculation;
    // 4b) Looping through the file
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
            if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants in " << durationOverall << "secs" << std::endl;
                std::cerr << "GettingCounts " << durationGettingCounts << " calculation " << durationCalculation << "secs" << std::endl;
            }
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1]; coordDouble = stringToDouble(coord);
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4]; bool ignoreSite = false;
            if (altAllele == "*") ignoreSite = true;
            if (!opt::allowIndels) {
                if (refAllele.length() > 1 || altAllele.length() > 1) ignoreSite = true;
            }
            if (ignoreSite) {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            startGettingCounts = std::clock();
            GeneralSetCounts* c = new GeneralSetCounts(popToPosMap, (int)genotypes.size());
            c->getSetVariantCountsSimple(genotypes, posToPopMap);
            genotypes.clear(); genotypes.shrink_to_fit();
            durationGettingCounts = ( std::clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            // std::cerr << "Here:" << totalVariantNumber << std::endl;
            
            startCalculation = std::clock();
            
           // for (std::vector<std::string>::size_type i = 0; i != populations.size(); i++) {
           //     allPs[i] = c->setAAFs.at(populations[i]);
           // }
            // durationFirstLoop = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            
            
            
            // Check if we are still in the same physical window...
            if (coordDouble > currentWindowEnd || coordDouble < currentWindowStart) {
                // Find the number of accessible BPs in this window:
                int accessibleInThisWindow = opt::fixedWindowSize;
                if (!opt::accesibleGenBedFile.empty()) {
                    if (coordDouble > currentWindowEnd) {
                        accessibleInThisWindow = ag->getAccessibleBPinRegion(chr, currentWindowStart, currentWindowStart + opt::fixedWindowSize);
                    } else if (coordDouble < currentWindowStart) {
                        accessibleInThisWindow = ag->getAccessibleBPinRegion(chr, 0, opt::fixedWindowSize);
                    }
                }
                
                for (int i = 0; i != outgroups.size(); i++) {
                    for (int j = 0; j != ingroups.size(); j++) {
                        int nSNPs = (int)DxyFixedWindowPerSNP[i][j].size();
                        if (nSNPs > 0) {
                            DxyFixedWindowAveraged[i][j] = vector_average_withRegion(DxyFixedWindowPerSNP[i][j], accessibleInThisWindow);
                        } else {
                            DxyFixedWindowAveraged[i][j] = 0.0;
                        }
                        DxyFixedWindowPerSNP[i][j].clear();
                    }
                    
                    *outFilesFixedWindow[i] << chr << "\t" << currentWindowStart << "\t" << currentWindowEnd << "\t" << usedVars[i] << "\t" << missingVars[i] << "\t" << accessibleInThisWindow << "\t";
                    print_vector(DxyFixedWindowAveraged[i], *outFilesFixedWindow[i]);
                    usedVars[i] = 0; missingVars[i] = 0;
                }
                
                if (coordDouble > currentWindowEnd) {
                    currentWindowStart = currentWindowStart + opt::fixedWindowSize; currentWindowEnd = currentWindowEnd + opt::fixedWindowSize;
                } else if (coordDouble < currentWindowStart) {
                    currentWindowStart = 0; currentWindowEnd = 0 + opt::fixedWindowSize;
                }
            }
            
            // std::cerr << coord << "\t";
            // print_vector_stream(SNPgeneDetails, std::cerr);
            
            // Now calculate Dxy per SNP:
            double AFout;
            for (int i = 0; i != outgroups.size(); i++) {
                AFout = c->setAAFs.at(outgroups[i]);
                if (AFout == -1) {
                    missingVars[i]++;
                    for (int j = 0; j != ingroups.size(); j++) { missingDist[i][j]++; }
                    continue;  // If the outgroup has entirely missing data, just move on to the next one
                }
                usedVars[i]++;
                for (int j = 0; j != ingroups.size(); j++) {
                    double AFin = c->setAAFs.at(ingroups[j]);
                    if (AFin != -1) {
                        double thisSNPdxy = DxyPerSNPfromAFs(AFout, c->setAAFs.at(ingroups[j]));
                        DxyFixedWindowPerSNP[i][j].push_back(thisSNPdxy);
                    } else {
                        missingDist[i][j]++;
                    }
                    
                }
            }
            durationCalculation = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    return 0;
    
    
    
    
    
}



void parseDistOutOptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'f': arg >> opt::fixedWindowSize; break;
            case 'n': arg >> opt::runName; break;
            case 'i': opt::allowIndels = true; break;
            case OPT_ACC_GEN_BED: arg >> opt::accesibleGenBedFile; break;
            case 'h':
                std::cout << DISTOUT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 4) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 4)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DISTOUT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::outgroupFile = argv[optind++];
    opt::ingroupFile = argv[optind++];
}
