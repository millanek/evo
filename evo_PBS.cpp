//
//  evo_PBS.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 03/01/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "evo_PBS.h"
#include "process_vcf_annotation_tools.h"
#include <deque>

#define SUBPROGRAM "PBS"

#define DEBUG 1

static const char *PBS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt PBS_trios.txt\n"
"Calculate the PBS statistic from:\n"
"Sequencing of 50 human exomes reveals adaptation to high altitude. Science 329, 75–78 (2010).\n"
"The POPULATIONS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The PBS_trios.txt should contain names of three populations for which the PBS will be calculated:\n"
"POP1   POP2    POP3\n"
"There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_POP3_PBS_SIZE_STEP.txt\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -w SIZE,STEP --window=SIZE,STEP         the parameters of the sliding window: contains SIZE SNPs and move by STEP (default: 20,10)\n"
"       --annot=ANNOTATION.gffExtract           (optional)gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs PBS per gene (only exons, with introns, and with 3kb upstream)\n"
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
    static string PBStriosFile;
    static string annotFile;
    static string runName = "";
    static int windowSize = 20;
    static int windowStep = 10;
}

inline std::vector<double> calculatePBSfromAFs(double p1, double p2, double p3, double p1AlleleCount, double p2AlleleCount, double p3AlleleCount) {
    // std::cerr << "p:\t" << p1 << "\t" << p2 << "\t" << p3 <<  std::endl;
    double Fst12; double Fst13; double Fst23;
    double power12 = pow(p1-p2, 2);
    double power13 = pow(p1-p3, 2);
    double power23 = pow(p2-p3, 2);
    double fraction1 = (p1*(1-p1))/(p1AlleleCount-1);
    double fraction2 = (p2*(1-p2))/(p2AlleleCount-1);
    double fraction3 = (p3*(1-p3))/(p3AlleleCount-1);
    double numerator12 = power12 - fraction1 - fraction2;
    double numerator13 = power13 - fraction1 - fraction3;
    double numerator23 = power23 - fraction2 - fraction3;
    double denominator12 = (p1*(1-p2))+(p2*(1-p1));
    double denominator13 = (p1*(1-p3))+(p3*(1-p1));
    double denominator23 = (p2*(1-p3))+(p3*(1-p2));
    if ((p1 == 0 && p2 == 0) || (p1 == 1 && p2 == 1)) { Fst12 = 0.0; } else { Fst12 = numerator12/denominator12; }
    if ((p1 == 0 && p3 == 0) || (p1 == 1 && p3 == 1)) { Fst13 = 0.0; } else { Fst13 = numerator13/denominator13; }
    if ((p2 == 0 && p3 == 0) || (p2 == 1 && p3 == 1)) { Fst23 = 0.0; } else { Fst23 = numerator23/denominator23; }
    // std::cerr << Fst12 << "\t" << Fst13 << "\t" << Fst23 <<  std::endl;
    if (Fst12 < 0) Fst12 = 0; if (Fst13 < 0) Fst13 = 0; if (Fst23 < 0) Fst23 = 0;
    if (Fst12 == 1) Fst12 = 1 - (Fst12/p1AlleleCount); if (Fst13 == 1) Fst13 = 1 - (Fst13/p1AlleleCount); if (Fst23 == 1) Fst23 = 1 - (Fst23/p2AlleleCount);
    double T12 = -log(1-Fst12); double T13 = -log(1-Fst13); double T23 = -log(1-Fst23);
    double PBS1 = (T12+T13-T23)/2.0;
    double PBS2 = (T12+T23-T13)/2.0;
    double PBS3 = (T13+T23-T12)/2.0;
    if (PBS1 < 0) PBS1 = 0; if (PBS2 < 0) PBS2 = 0; if (PBS3 < 0) PBS3 = 0;
    //std::cerr << PBS1 << "\t" << PBS2 << "\t" << PBS3 <<  std::endl;
    std::vector<double> PBS; PBS.push_back(PBS1); PBS.push_back(PBS2); PBS.push_back(PBS3);
    return PBS;
}


inline std::vector<double> calculatePBSnumerators(const GeneralSetCounts* c,std::vector<string>& PBStrio) {
    double power12 = pow((c->setAAFs.at(PBStrio[0])-c->setAAFs.at(PBStrio[1])), 2);
    double power13 = pow((c->setAAFs.at(PBStrio[0])-c->setAAFs.at(PBStrio[2])), 2);
    double power23 = pow((c->setAAFs.at(PBStrio[1])-c->setAAFs.at(PBStrio[2])), 2);
    double fraction1 = (c->setAAFs.at(PBStrio[0])*(1-c->setAAFs.at(PBStrio[0])))/(c->setAlleleCounts.at(PBStrio[0])-1);
    double fraction2 = (c->setAAFs.at(PBStrio[1])*(1-c->setAAFs.at(PBStrio[1])))/(c->setAlleleCounts.at(PBStrio[1])-1);
    double fraction3 = (c->setAAFs.at(PBStrio[2])*(1-c->setAAFs.at(PBStrio[2])))/(c->setAlleleCounts.at(PBStrio[2])-1);
    double numerator12 = power12 - fraction1 - fraction2;
    double numerator13 = power13 - fraction1 - fraction3;
    double numerator23 = power23 - fraction2 - fraction3;
    std::vector<double> PBSnumerators; PBSnumerators.push_back(numerator12);
    PBSnumerators.push_back(numerator13); PBSnumerators.push_back(numerator23);
    return PBSnumerators;
}

inline std::vector<double> calculatePBSdenominator(const GeneralSetCounts* c,std::vector<string>& PBStrio) {
    double denominator12 = (c->setAAFs.at(PBStrio[0])*(1-c->setAAFs.at(PBStrio[1])))+(c->setAAFs.at(PBStrio[1])*(1-c->setAAFs.at(PBStrio[0])));
    double denominator13 = (c->setAAFs.at(PBStrio[0])*(1-c->setAAFs.at(PBStrio[2])))+(c->setAAFs.at(PBStrio[2])*(1-c->setAAFs.at(PBStrio[0])));
    double denominator23 = (c->setAAFs.at(PBStrio[1])*(1-c->setAAFs.at(PBStrio[2])))+(c->setAAFs.at(PBStrio[2])*(1-c->setAAFs.at(PBStrio[1])));
    std::vector<double> PBSdenominators; PBSdenominators.push_back(denominator12);
    PBSdenominators.push_back(denominator13); PBSdenominators.push_back(denominator23);
    return PBSdenominators;
}


int PBSmain(int argc, char** argv) {
    parsePBSoptions(argc, argv);
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
    std::ifstream* PBStriosFile = new std::ifstream(opt::PBStriosFile.c_str());
    
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
    
    // Get the PBS trios
    std::vector<std::ofstream*> outFiles;
    std::vector<std::ofstream*> outFilesGenes;
    std::vector<std::vector<string> > PBStrios;
    while (getline(*PBStriosFile,line)) {
        // std::cerr << line << std::endl;
        std::vector<string> threePops = split(line, '\t'); assert(threePops.size() == 3);
        std::ofstream* outFile = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBS_" + opt::runName + "_" + numToString(opt::windowSize) + "_" + numToString(opt::windowStep) + ".txt");
        *outFile << "chr\t" << threePops[0] << "\t" << threePops[1] << "\t" << threePops[2] << std::endl;
        outFiles.push_back(outFile);
        if (!opt::annotFile.empty()) {
            std::ofstream* outFileGenes = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_PBSGenes_" + opt::runName + "_" + numToString(opt::windowSize) + "_" + numToString(opt::windowStep) + ".txt");
           // *outFileGenes << "gene\t" << "numSNPsExons\t" << "numSNPsWithIntrons\t" << "numSNPsWith3kbUpstr\t" << threePops[0] << "_exons\t" << threePops[1] << "_exons\t" << threePops[2] << "_exons\t" << threePops[0] << "_wIntrons\t" << threePops[1] << "_wIntrons\t" << threePops[2] << "_wIntrons\t" << threePops[0] << "_w3kbUpstr\t" << threePops[1] << "_w3kbUpstr\t" << threePops[2] << "_w3kbUpstr" << std::endl;
            *outFileGenes << "gene\t" << "numSNPsExons\t" << "numSNPsWithIntrons\t" << threePops[0] << "_exons\t" << threePops[1] << "_exons\t" << threePops[2] << "_exons\t" << threePops[0] << "_wIntrons\t" << threePops[1] << "_wIntrons\t" << threePops[2] << "_wIntrons" << std::endl;
            outFilesGenes.push_back(outFileGenes);
        }
        PBStrios.push_back(threePops);
    }
    
    // And need to prepare the vectors to hold the PBS values and the coordinates:
    std::deque<double> initDeq(opt::windowSize,0.0); // deque to initialise per-site PBS values
    std::vector<std::deque<double>> initThreeDeques(4,initDeq); // vector of three per-site PBS deques - for each population in the trio, and the fourth is for the coordinates
    std::vector<std::vector<std::deque<double>>> PBSresults(PBStrios.size(),initThreeDeques);
    
    // if (!opt::annotFile.empty()) {
    std::vector<std::vector<double>> initGeneVectors(9); // For the nine PBS columns in the _PBSGenes_ files
    std::vector<std::vector<std::vector<double>>> PBSgeneResults(PBStrios.size(),initGeneVectors);
    std::string currentGene = "";
    //}
    //std::deque<std::vector<double>> regionPBSnums; regionPBSnums.assign(opt::windowSize,init);
    //std::deque<std::vector<double>> regionPBSdenoms; regionPBSdenoms.assign(opt::windowSize,init);
    //std::vector<double> allPs(populations.size(),0.0);
    int totalVariantNumber = 0;
    std::vector<int> usedVars(PBStrios.size(),0); // Will count the number of used variants for each trio
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
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1) {
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
            double p1; double p2; double p3;
            for (int i = 0; i != PBStrios.size(); i++) {
                p1 = c->setAAFs.at(PBStrios[i][0]); //assert(p_S1 == pS1test);
                if (p1 == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                p2 = c->setAAFs.at(PBStrios[i][1]); //assert(p_S2 == pS2test);
                if (p2 == -1) continue;
                p3 = c->setAAFs.at(PBStrios[i][2]); // assert(p_S3 == pS3test);
                if (p3 == -1) continue;
                
                if (p1 == 0 && p2 == 0 && p3 == 0) { continue; }
                if (p1 == 1 && p2 == 1 && p3 == 1) { continue; }
                usedVars[i]++;
                
                std::vector<double> thisSNP_PBS = calculatePBSfromAFs(p1,p2,p3,
                                                                      c->setAlleleCounts.at(PBStrios[i][0]),
                                                                      c->setAlleleCounts.at(PBStrios[i][1]),
                                                                      c->setAlleleCounts.at(PBStrios[i][2]));
                
                PBSresults[i][0].push_back(thisSNP_PBS[0]); PBSresults[i][1].push_back(thisSNP_PBS[1]); PBSresults[i][2].push_back(thisSNP_PBS[2]);
                PBSresults[i][3].push_back(stringToDouble(coord));
                PBSresults[i][0].pop_front(); PBSresults[i][1].pop_front(); PBSresults[i][2].pop_front(); PBSresults[i][3].pop_front();
                
                if (!opt::annotFile.empty()) { if (SNPgeneDetails[0] != "") {
                    if (SNPgeneDetails[1] == "exon") {
                        PBSgeneResults[i][0].push_back(thisSNP_PBS[0]); PBSgeneResults[i][1].push_back(thisSNP_PBS[1]); PBSgeneResults[i][2].push_back(thisSNP_PBS[2]);
                    } PBSgeneResults[i][3].push_back(thisSNP_PBS[0]); PBSgeneResults[i][4].push_back(thisSNP_PBS[1]); PBSgeneResults[i][5].push_back(thisSNP_PBS[2]);
                    currentGene = SNPgeneDetails[0];
                }}
                if (usedVars[i] > opt::windowSize && (usedVars[i] % opt::windowStep == 0)) {
                    // std::cerr << PBSresults[i][0][0] << std::endl;
                    *outFiles[i] << chr << "\t" << PBSresults[i][3][0] << "\t" << coord << "\t" << vector_average(PBSresults[i][0]) << "\t" << vector_average(PBSresults[i][1]) << "\t" << vector_average(PBSresults[i][2]) << std::endl;
                }
                // }
            }
            if (!opt::annotFile.empty()) { if (SNPgeneDetails[0] == "" && currentGene != "") {
                for (int i = 0; i != PBStrios.size(); i++) {
                    *outFilesGenes[i] << currentGene << "\t" << PBSgeneResults[i][0].size() << "\t" << PBSgeneResults[i][3].size() << "\t" << vector_average(PBSgeneResults[i][0]) << "\t" << vector_average(PBSgeneResults[i][1]) << "\t" << vector_average(PBSgeneResults[i][2]) << "\t" << vector_average(PBSgeneResults[i][3]) << "\t" << vector_average(PBSgeneResults[i][4]) << "\t" << vector_average(PBSgeneResults[i][5]) << std::endl;
                        for (int j = 0; j <= 5; j++) { PBSgeneResults[i][j].clear(); }
                }
                currentGene = "";
            }}
            durationCalculation = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    return 0;
    
}



void parsePBSoptions(int argc, char** argv) {
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
                std::cout << PBS_USAGE_MESSAGE;
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
        std::cout << "\n" << PBS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::PBStriosFile = argv[optind++];
}
