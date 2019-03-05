//
//  evo_ABS.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 01/03/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "evo_ABS.h"
#include "process_vcf_annotation_tools.h"
#include <deque>

#define SUBPROGRAM "ABS"

#define DEBUG 1

static const char *ABS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt ABS_quartets.txt\n"
"Calculate the ABS statistic from:\n"
"Cheng, X., Xu, C. & Degiorgio, M. Fast and robust detection of ancestral selective sweeps. Mol Ecol 26, 6871–6891 (2017).\n"
"The POPULATIONS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The ABS_quartets.txt should names of the four populations for which the ABS will be calculated, in the following order:\n"
"POP_W   POP_X    POP_Y POP_Z\n"
"There can be multiple lines. The program generates one ouput file, named like RUN_NAME_ABS_SIZE_STEP.txt\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -w SIZE,STEP --window=SIZE,STEP         the parameters of the sliding window: contains SIZE SNPs and move by STEP (default: 20,10)\n"
"       --annot=ANNOTATION.gffExtract           (optional)gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs ABS per gene (only exons, with introns, and with 3kb upstream)\n"
"       -r , --region=start,length              (optional) only process a subset of the VCF file\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hw:n:";

enum { OPT_ANNOT  };

static const struct option longopts[] = {
    { "window",   required_argument, NULL, 'w' },
    { "annot",   required_argument, NULL, OPT_ANNOT },
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string ABSquartetFile;
    static string annotFile;
    static string runName = "";
    static int windowSize = 20;
    static int windowStep = 10;
}

inline std::vector<double> calculateABSfromAFs(double pW, double pX, double pY, double pZ,
                                               double pWAlleleCount, double pXAlleleCount, double pYAlleleCount, double pZAlleleCount) {
    // std::cerr << "p:\t" << p1 << "\t" << p2 << "\t" << p3 <<  std::endl;
    double FstWX; double FstWY; double FstWZ;
    double FstXY; double FstXZ; double FstYZ;
    double powerWX = pow(pW-pX, 2); double powerWY = pow(pW-pY, 2);
    double powerWZ = pow(pW-pZ, 2); double powerXY = pow(pX-pY, 2);
    double powerXZ = pow(pX-pZ, 2); double powerYZ = pow(pY-pZ, 2);
    double fractionX = (pW*(1-pW))/(pWAlleleCount-1); double fractionW = (pX*(1-pX))/(pXAlleleCount-1);
    double fractionY = (pY*(1-pY))/(pYAlleleCount-1); double fractionZ = (pZ*(1-pZ))/(pZAlleleCount-1);
    double numeratorWX = powerWX - fractionW - fractionX;
    double numeratorWY = powerWY - fractionW - fractionY;
    double numeratorWZ = powerWZ - fractionW - fractionZ;
    double numeratorXY = powerXY - fractionX - fractionY;
    double numeratorXZ = powerXZ - fractionX - fractionZ;
    double numeratorYZ = powerYZ - fractionY - fractionZ;
    double denominatorWX = (pW*(1-pX))+(pX*(1-pW)); double denominatorWY = (pW*(1-pY))+(pY*(1-pW));
    double denominatorWZ = (pW*(1-pZ))+(pZ*(1-pW)); double denominatorXY = (pX*(1-pY))+(pY*(1-pX));
    double denominatorXZ = (pX*(1-pZ))+(pZ*(1-pX)); double denominatorYZ = (pY*(1-pZ))+(pZ*(1-pY));
    if ((pW == 0 && pX == 0) || (pW == 1 && pX == 1)) { FstWX = 0.0; } else { FstWX = numeratorWX/denominatorWX; }
    if ((pW == 0 && pY == 0) || (pW == 1 && pY == 1)) { FstWY = 0.0; } else { FstWY = numeratorWY/denominatorWY; }
    if ((pW == 0 && pZ == 0) || (pW == 1 && pZ == 1)) { FstWZ = 0.0; } else { FstWZ = numeratorWZ/denominatorWZ; }
    if ((pX == 0 && pY == 0) || (pX == 1 && pY == 1)) { FstXY = 0.0; } else { FstXY = numeratorXY/denominatorXY; }
    if ((pX == 0 && pZ == 0) || (pX == 1 && pZ == 1)) { FstXZ = 0.0; } else { FstXZ = numeratorXZ/denominatorXZ; }
    if ((pY == 0 && pZ == 0) || (pY == 1 && pZ == 1)) { FstYZ = 0.0; } else { FstYZ = numeratorYZ/denominatorYZ; }
    // std::cerr << Fst12 << "\t" << Fst13 << "\t" << Fst23 <<  std::endl;
    if (FstWX < 0) FstWX = 0; if (FstWY < 0) FstWY = 0; if (FstWZ < 0) FstWZ = 0;
    if (FstXY < 0) FstXY = 0; if (FstXZ < 0) FstXZ = 0; if (FstYZ < 0) FstYZ = 0;
    if (FstWX == 1) FstWX = 1 - (FstWX/((pWAlleleCount+pXAlleleCount)/2.0));
    if (FstWY == 1) FstWY = 1 - (FstWY/((pWAlleleCount+pYAlleleCount)/2.0));
    if (FstWZ == 1) FstWZ = 1 - (FstWZ/((pWAlleleCount+pZAlleleCount)/2.0));
    if (FstXY == 1) FstXY = 1 - (FstXY/((pXAlleleCount+pYAlleleCount)/2.0));
    if (FstXZ == 1) FstXZ = 1 - (FstXZ/((pXAlleleCount+pZAlleleCount)/2.0));
    if (FstYZ == 1) FstYZ = 1 - (FstYZ/((pYAlleleCount+pZAlleleCount)/2.0));
    double TWX = -log(1-FstWX); double TWY = -log(1-FstWY); double TWZ = -log(1-FstWZ);
    double TXY = -log(1-FstXY); double TXZ = -log(1-FstXZ); double TYZ = -log(1-FstYZ);
    double TWYTXZ = TWY + TXZ; double TWZTXY = TWZ + TXY; double TWXTYZ = TWX + TYZ;
    double twoMax; if (TWYTXZ >= TWZTXY) { twoMax = TWYTXZ; } else { twoMax = TWZTXY; }
    double ABSalt = (twoMax-TWXTYZ)/2.0;
    double threeMax; double threeMin;
    if (TWYTXZ >= TWZTXY && TWYTXZ >= TWXTYZ) { threeMax = TWYTXZ; if (TWZTXY <= TWXTYZ) { threeMin = TWZTXY;} else { threeMin = TWXTYZ;}}
    else if (TWZTXY >= TWYTXZ && TWZTXY >= TWXTYZ) { threeMax = TWZTXY; if (TWYTXZ <= TWXTYZ) { threeMin = TWYTXZ;} else { threeMin = TWXTYZ;}}
    else { threeMax = TWXTYZ; if (TWYTXZ <= TWZTXY) { threeMin = TWYTXZ;} else { threeMin = TWZTXY;} }
    double ABSmain = (threeMax-threeMin)/2.0;
    //if (PBS1 < 0) PBS1 = 0; if (PBS2 < 0) PBS2 = 0; if (PBS3 < 0) PBS3 = 0;
    //std::cerr << PBS1 << "\t" << PBS2 << "\t" << PBS3 <<  std::endl;
    std::vector<double> ABS; ABS.push_back(ABSmain); ABS.push_back(ABSalt);
    return ABS;
}


int ABSmain(int argc, char** argv) {
    parseABSoptions(argc, argv);
    string line; // for reading the input files
    
    // Load up the annotation file
    Annotation wgAnnotation; std::ifstream* annotFile;
    if (!opt::annotFile.empty()) {
        annotFile = new std::ifstream(opt::annotFile.c_str());
        Annotation Annot(annotFile, false); // Does not use transcripts annotated as 5' or 3' partial
        wgAnnotation = Annot;
    }
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    std::ifstream* ABSquartetFile = new std::ifstream(opt::ABSquartetFile.c_str());
    
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
    
    // Get the ABS quartets
    std::vector<std::ofstream*> outFiles;
    //std::vector<std::ofstream*> outFileGenes;
    std::vector<std::vector<string> > ABSquartets;
    while (getline(*ABSquartetFile,line)) {
        // std::cerr << line << std::endl;
        std::vector<string> fourPops = split(line, '\t'); assert(fourPops.size() == 4);
        ABSquartets.push_back(fourPops);
        std::ofstream* outFile = new std::ofstream(fourPops[0] + "_" + fourPops[1] + "_" + fourPops[2] + "_" + fourPops[3] + "_ABS_" + opt::runName + "_" + numToString(opt::windowSize) + "_" + numToString(opt::windowStep) + ".txt");
        *outFile << "chr\tpos1\tpos2\tABS\tABSalt" << std::endl;
        outFiles.push_back(outFile);
    }
    
    // And need to prepare the vectors to hold the ABS values and the coordinates:
    std::deque<double> initDeq(opt::windowSize,0.0); // deque to initialise per-site ABS values
    std::vector<std::deque<double>> initABSdeques(3,initDeq); // vector of per-site ABS deques - two for each quartet plus one for the coordinates
    std::vector<std::vector<std::deque<double>>> ABSresults(ABSquartets.size(),initABSdeques); // one set of deques for each quartet
    
    // if (!opt::annotFile.empty()) {
    //std::vector<std::vector<double>> initGeneVectors(9); // For the nine PBS columns in the _PBSGenes_ files
    //std::vector<std::vector<std::vector<double>>> PBSgeneResults(PBStrios.size(),initGeneVectors);
    //std::string currentGene = "";
    //}
    //std::deque<std::vector<double>> regionPBSnums; regionPBSnums.assign(opt::windowSize,init);
    //std::deque<std::vector<double>> regionPBSdenoms; regionPBSdenoms.assign(opt::windowSize,init);
    //std::vector<double> allPs(populations.size(),0.0);
    int totalVariantNumber = 0;
    std::vector<int> usedVars(ABSquartets.size(),0); // Will count the number of used variants for each quartet
    std::vector<string> sampleNames; std::vector<std::string> fields;
    int reportProgressEvery = 1000; string chr; string coord;
    std::clock_t start; std::clock_t startGettingCounts; std::clock_t startCalculation;
    double durationOverall; double durationGettingCounts; double durationCalculation;
    
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
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1];
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            //std::vector<std::string> info = split(fields[7], ';');
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1 || altAllele == "*") {
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
            
            // find if we are in a gene:
            std::vector<string> SNPgeneDetails = wgAnnotation.getSNPgeneDetails(chr, atoi(coord.c_str()));
            
            // std::cerr << coord << "\t";
            // print_vector_stream(SNPgeneDetails, std::cerr);
            // Now calculate the PBS stats:
            double pW; double pX; double pY; double pZ;
            for (int i = 0; i != ABSquartets.size(); i++) {
                pW = c->setAAFs.at(ABSquartets[i][0]);
                if (pW == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                pX = c->setAAFs.at(ABSquartets[i][1]);
                if (pX == -1) continue;
                pY = c->setAAFs.at(ABSquartets[i][2]);
                if (pY == -1) continue;
                pZ = c->setAAFs.at(ABSquartets[i][3]);
                if (pZ == -1) continue;
                
                if (pW == 0 && pX == 0 && pY == 0 && pZ == 0) { continue; }
                if (pW == 1 && pX == 1 && pY == 1 && pZ == 1) { continue; }
                usedVars[i]++;
                
                std::vector<double> thisSNP_ABS = calculateABSfromAFs(pW,pX,pY,pZ,
                                                                      c->setAlleleCounts.at(ABSquartets[i][0]),
                                                                      c->setAlleleCounts.at(ABSquartets[i][1]),
                                                                      c->setAlleleCounts.at(ABSquartets[i][2]),
                                                                      c->setAlleleCounts.at(ABSquartets[i][3]));
                
                ABSresults[i][0].push_back(thisSNP_ABS[0]); ABSresults[i][1].push_back(thisSNP_ABS[1]); ABSresults[i][2].push_back(stringToDouble(coord));
                ABSresults[i][0].pop_front(); ABSresults[i][1].pop_front(); ABSresults[i][2].pop_front();
                
                /*if (!opt::annotFile.empty()) { if (SNPgeneDetails[0] != "") {
                    if (SNPgeneDetails[1] == "exon") {
                        PBSgeneResults[i][0].push_back(thisSNP_PBS[0]); PBSgeneResults[i][1].push_back(thisSNP_PBS[1]); PBSgeneResults[i][2].push_back(thisSNP_PBS[2]);
                    } PBSgeneResults[i][3].push_back(thisSNP_PBS[0]); PBSgeneResults[i][4].push_back(thisSNP_PBS[1]); PBSgeneResults[i][5].push_back(thisSNP_PBS[2]);
                    currentGene = SNPgeneDetails[0];
                }} */
                if (usedVars[i] > opt::windowSize && (usedVars[i] % opt::windowStep == 0)) {
                    // std::cerr << PBSresults[i][0][0] << std::endl;
                    *outFiles[i] << chr << "\t" << ABSresults[i][2][0] << "\t" << coord << "\t" << vector_average(ABSresults[i][0]) << "\t" << vector_average(ABSresults[i][1]) << std::endl;
                }
                // }
            }
            /*if (!opt::annotFile.empty()) { if (SNPgeneDetails[0] == "" && currentGene != "") {
                for (int i = 0; i != ABStrios.size(); i++) {
                    *outFilesGenes[i] << currentGene << "\t" << ABSgeneResults[i][0].size() << "\t" << PBSgeneResults[i][3].size() << "\t" << vector_average(PBSgeneResults[i][0]) << "\t" << vector_average(PBSgeneResults[i][1]) << "\t" << vector_average(PBSgeneResults[i][2]) << "\t" << vector_average(PBSgeneResults[i][3]) << "\t" << vector_average(PBSgeneResults[i][4]) << "\t" << vector_average(PBSgeneResults[i][5]) << std::endl;
                    for (int j = 0; j <= 5; j++) { PBSgeneResults[i][j].clear(); }
                }
                currentGene = "";
            }} */
            durationCalculation = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    return 0;
    
}



void parseABSoptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case 'n': arg >> opt::runName; break;
            case OPT_ANNOT: arg >> opt::annotFile; break;
            case 'h':
                std::cout << ABS_USAGE_MESSAGE;
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
        std::cout << "\n" << ABS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::ABSquartetFile = argv[optind++];
}

