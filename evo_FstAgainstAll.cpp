//
//  evo_FstAgainstAll.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 07/02/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "evo_FstAgainstAll.h"


#include "process_vcf_annotation_tools.h"
#include "process_vcf_fst.h"
#include <deque>

#define SUBPROGRAM "FstGlobal"

#define DEBUG 1

static const char *FSTGLOBAL_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf POPULATIONS.txt\n"
"Calculate the Fst statistic for each population from the POPULATIONS.txt file, when compared against the rest of the dataset:\n"
"The POPULATIONS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"samples under POPULATION_IDs that are either \"Outgroup\" or \"xxx\" will not be included in the 'rest of the dataset'\n"
"The program generates one ouput file, named like RUN_NAME_FstGlobal_SIZE_STEP.txt\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -w SIZE,STEP --window=SIZE,STEP         the parameters of the sliding window: contains SIZE SNPs and move by STEP (default: 20,10)\n"
"       --annot=ANNOTATION.gffExtract           (optional)gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs PBS per gene (only exons, with introns, and with 3kb upstream)\n"
"       -r , --region=start,length              (optional) only process a subset of the VCF file\n"
"       -n, --run-name=RUN_NAME                 run-name will be included in the output file name\n"
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

int FstGlobalMain(int argc, char** argv) {
    parseFstGlobaloptions(argc, argv);
    string line; // for reading the input files
    
    // Load up the annotation file
    Annotation wgAnnotation; std::ifstream* annotFile;
    if (!opt::annotFile.empty()) {
        annotFile = new std::ifstream(opt::annotFile.c_str());
        Annotation Annot(annotFile, false); // Does not use transcripts annotated as 5' or 3' partial
        wgAnnotation = Annot;
    }
    std::vector<std::vector<string> > annotation;
    
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
    std::vector<string> populations; std::vector<string> populationsToUse;
    for(std::map<string,std::vector<string>>::iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
        populations.push_back(it->first);
        if (it->first != "Outgroup" && it->first != "xxx") {
            populationsToUse.push_back(it->first);
        }
        // std::cerr << it->first << std::endl;
    } std::cerr << "There are " << populations.size() << " populations " << std::endl;
    
    std::ofstream* outFile = new std::ofstream(opt::runName + "_FstGlobal_" + numToString(opt::windowSize) + "_" + numToString(opt::windowStep) + ".txt");
    *outFile << "chr\tpos1\tpos2\tFstGlobal"; for (int i = 0; i < populationsToUse.size(); i++) { *outFile << "\t" << populationsToUse[i];}
    *outFile << std::endl;
    if (!opt::annotFile.empty()) {
        std::ofstream* outFileGenes = new std::ofstream(opt::runName + "_FstGlobalGenes_" + opt::runName + "_" + numToString(opt::windowSize) + "_" + numToString(opt::windowStep) + ".txt");
        // *outFileGenes << "gene\t" << "numSNPsExons\t" << "numSNPsWithIntrons\t" << "numSNPsWith3kbUpstr\t" << threePops[0] << "_exons\t" << threePops[1] << "_exons\t" << threePops[2] << "_exons\t" << threePops[0] << "_wIntrons\t" << threePops[1] << "_wIntrons\t" << threePops[2] << "_wIntrons\t" << threePops[0] << "_w3kbUpstr\t" << threePops[1] << "_w3kbUpstr\t" << threePops[2] << "_w3kbUpstr" << std::endl;
        *outFileGenes << "gene\t" << "numSNPsExons\t" << "numSNPsWithIntrons";
        for (int i = 0; i < populationsToUse.size(); i++) {
            *outFileGenes << "\t" << populationsToUse[i] << "_exons";
            *outFileGenes << "\t" << populationsToUse[i] << "_wIntrons";
        } *outFileGenes << std::endl;
    }
    
    // And need to prepare the vectors to hold the Fst values and the coordinates:
    std::deque<double> initFst(opt::windowSize,0.0); // deque to initialise per-site PBS values
    std::vector<std::deque<double>> FstNumDeques(populationsToUse.size()+1,initFst); // Fst numerators for each population, and one for the global Fst
    std::vector<std::deque<double>> FstDenomDeques(populationsToUse.size()+1,initFst); // Fst denominators for each population, and one for the global Fst
    std::deque<string> coordDeque(opt::windowSize,"0");
    
    // if (!opt::annotFile.empty()) {
    std::vector<std::vector<double>> FstGeneNumVectors((2*populationsToUse.size())+2); // For the nine PBS columns in the _PBSGenes_ files
    std::vector<std::vector<double>> FstGeneDenumVectors((2*populationsToUse.size())+2); // For the nine PBS columns in the _PBSGenes_ files
    std::string currentGene = "NN";
    //}
    //std::deque<std::vector<double>> regionPBSnums; regionPBSnums.assign(opt::windowSize,init);
    //std::deque<std::vector<double>> regionPBSdenoms; regionPBSdenoms.assign(opt::windowSize,init);
    //std::vector<double> allPs(populations.size(),0.0);
    int totalVariantNumber = 0; int usedVariantNumber = 0;
    //std::vector<int> usedVars(PBStrios.size(),0); // Will count the number of used variants for each trio
    std::vector<string> sampleNames; std::vector<std::string> fields;
    int reportProgressEvery = 10000; string chr; string coord;
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
            usedVariantNumber++;
            startGettingCounts = std::clock();
            GeneralSetCountsWithComplements* c = new GeneralSetCountsWithComplements(popToPosMap, (int)genotypes.size());
            c->getSetVariantCountsSimple(genotypes, posToPopMap);
            c->getComplementCounts(populationsToUse);
            genotypes.clear(); genotypes.shrink_to_fit();
            durationGettingCounts = ( std::clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            // std::cerr << "Here:" << totalVariantNumber << std::endl;
            
            startCalculation = std::clock();
            
            
            // for (std::vector<std::string>::size_type i = 0; i != populations.size(); i++) {
            //     allPs[i] = c->setAAFs.at(populations[i]);
            // }
            // durationFirstLoop = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            
            // find if we are in a gene:
            std::vector<string> SNPgeneDetails = wgAnnotation.getSNPgeneDetails(chr, coord);
            
            coordDeque.push_back(coord); coordDeque.pop_front();
            *outFile << chr << "\t" << coordDeque[0] << "\t" << coord << "\t" << "N"; double Fst;
            // Now calculate the Fst stats:
            double p1; double p2; int n1; int n2;
            for (int i = 0; i != populationsToUse.size(); i++) {
                p1 = c->setAAFs.at(populationsToUse[i]); //assert(p_S1 == pS1test);
                p2 = c->setAAFsComplement.at(populationsToUse[i]);
                n1 = c->setAlleleCounts.at(populationsToUse[i]);
                n2 = c->setAlleleCountsComplement.at(populationsToUse[i]);
                // std::cerr << p1 << "\t" << p2 << "\t" << n1 << "\t" << n2 << std::endl;
                
                FstNumDeques[i].push_back(calculateFstNumerator(p1, p2, n1, n2)); FstNumDeques[i].pop_front();
                FstDenomDeques[i].push_back(calculateFstDenominator(p1, p2)); FstDenomDeques[i].pop_front();
                
               /* if (!opt::annotFile.empty()) { if (SNPgeneDetails[0] != "") {
                    if (SNPgeneDetails[0] == currentGene) {
                        if (SNPgeneDetails[1] == "exon") {
                            PBSgeneResults[i][0].push_back(thisSNP_PBS[0]); PBSgeneResults[i][1].push_back(thisSNP_PBS[1]); PBSgeneResults[i][2].push_back(thisSNP_PBS[2]);
                        } PBSgeneResults[i][3].push_back(thisSNP_PBS[0]); PBSgeneResults[i][4].push_back(thisSNP_PBS[1]); PBSgeneResults[i][5].push_back(thisSNP_PBS[2]);
                    } else {
                        if (currentGene != "NN") {
                            *outFilesGenes[i] << currentGene << "\t" << PBSgeneResults[i][0].size() << "\t" << PBSgeneResults[i][3].size() << "\t" << vector_average(PBSgeneResults[i][0]) << "\t" << vector_average(PBSgeneResults[i][1]) << "\t" << vector_average(PBSgeneResults[i][2]) << "\t" << vector_average(PBSgeneResults[i][3]) << "\t" << vector_average(PBSgeneResults[i][4]) << "\t" << vector_average(PBSgeneResults[i][5]) << std::endl;
                            for (int j = 0; j <= 5; j++) { PBSgeneResults[i][j].clear(); }
                        } currentGene = SNPgeneDetails[0];
                    }
                }} */
                if ((usedVariantNumber > opt::windowSize || opt::windowSize == opt::windowStep) && (usedVariantNumber % opt::windowStep == 0)) {
                    // std::cerr << PBSresults[i][0][0] << std::endl;
                    Fst = vector_average(FstNumDeques[i])/vector_average(FstDenomDeques[i]); if (Fst < 0 || vector_average(FstDenomDeques[i]) == 0) { Fst = 0; }
                    *outFile << "\t" << Fst;
                }
                // }
            }
            *outFile << std::endl;
            durationCalculation = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    return 0;
    
}



void parseFstGlobaloptions(int argc, char** argv) {
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
                std::cout << FSTGLOBAL_USAGE_MESSAGE;
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
        std::cout << "\n" << FSTGLOBAL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}
