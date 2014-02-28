//
//  process_vcf_cbs.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 22/11/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_cbs.h"
/*
 TO DO:
 High priority:
 - deal with the possibility that a scaffold may not have any variants (DONE - apart from scaffold_0)
 
 */


#define SUBPROGRAM "cbs"

#define DEBUG 1

static const char *CBS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE(.vcf or .fa, depending on options)\n"
"Calculate the lenghts of tracts of 'Compatibility by sequence' (cbs) between samples from genotypes in a VCF file\n"
"Also does some related tasks\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       --prepare-genome                        find the coordinates of undetermined bases in the reference genome\n"
"                                               if this option is used, INPUT_FILE.fa should be a fasta file with the reference genome\n"
"                                               and output is written into INPUT_FILE.ns (zero indexed coordinates)\n"
"       --cbs=UNDETERMIED.ns                    Calculate the lenghts of cbs tracts, taking into account the location of undetermined characters\n"
"                                               (Ns) in the reference genome, supplied in the UNDETERMIED.ns file (zero indexed coordinates)\n"
"       -m MIN, --min-sc-length=MIN             Only analyse genomic scaffolds of length >= MIN (default: MIN=0)\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the output\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_CBS, OPT_PREPARE_GENOME };

static const char* shortopts = "hs:m:";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "min-sc-length",   required_argument, NULL, 'm' },
    { "cbs",   required_argument, NULL, OPT_CBS },
    { "prepare-genome",   no_argument, NULL, OPT_PREPARE_GENOME },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genomeFile;
    static string undeterminedLociFile;
    static bool bPrepareGenome = false;
    static string sampleNameFile;
    static int minScLength = 0;
}



void findUndeterminedRegions() {
    // Open a connection to read from the genome file
    std::ifstream* genomeFile = new std::ifstream(opt::genomeFile.c_str());
    string genomeFileRoot = stripExtension(opt::genomeFile);
    string undeterminedRegionsFileName = genomeFileRoot + ".ns";
    std::ofstream* undeterminedRegionsFile = new std::ofstream(undeterminedRegionsFileName.c_str());
    
    string nextScaffoldName;
    getline(*genomeFile, nextScaffoldName);
    string currentScaffoldReference;
    int thisUndeterminedRegionStart;
    int thisUndeterminedRegionEnd;
    while (nextScaffoldName != "") {
        std::cerr << "Processing " << nextScaffoldName << "..." << std::endl;
        string thisScaffold = nextScaffoldName;
        currentScaffoldReference = readScaffold(genomeFile, nextScaffoldName);
        *undeterminedRegionsFile << thisScaffold << "\t" << currentScaffoldReference.length() << std::endl;
        bool bInUndetermined = false;
        for (string::size_type i = 0; i != currentScaffoldReference.size(); i++) {
            if (currentScaffoldReference[i] == 'N') {
                if (!bInUndetermined)
                    thisUndeterminedRegionStart = (int)i;
                if (i == currentScaffoldReference.size()-1) {
                    if (bInUndetermined)
                        *undeterminedRegionsFile << thisUndeterminedRegionStart << "\t" << i << std::endl;
                    else
                        *undeterminedRegionsFile << i << "\t" << i << std::endl;
                }
                bInUndetermined = true;
            } else {
                if (bInUndetermined) {
                    thisUndeterminedRegionEnd = (int)i-1;
                    *undeterminedRegionsFile << thisUndeterminedRegionStart << "\t" << thisUndeterminedRegionEnd << std::endl;
                    bInUndetermined = false;
                }
            }
        }
    }
    
    string line;
}

std::vector<string> numToStringVector(std::vector<int>& vec) {
    std::vector<string> converted;
    for (std::vector<int>::size_type i = 0; i != vec.size(); i++) {
        converted.push_back(numToString(vec[i]));
    }
    return converted;
}

std::vector<string> numToStringMinFilterVector(std::vector<int>& vec, int min) {
    std::vector<string> converted;
    for (std::vector<int>::size_type i = 0; i != vec.size(); i++) {
        if (vec[i] >= min)
            converted.push_back(numToString(vec[i]));
    }
    return converted;
}




bool compareVectorsByLength (std::vector<string> a,std::vector<string> b) {
    if (a.size() >= b.size()) {
        return true;
    } else {
        return false;
    }
}

void calculateCbs() {
    // Open a connection to read from the vcf file
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* unlFile = new std::ifstream(opt::undeterminedLociFile.c_str());
    string vcfFileRoot = stripExtension(opt::vcfFile);
    string cbsFileName = vcfFileRoot + ".cbsTracts";
    std::ofstream* cbsFile = new std::ofstream(cbsFileName.c_str());
    string cbsM10000FileName = vcfFileRoot + ".cbsTractsMin10000";
    std::ofstream* cbsFileM10000 = new std::ofstream(cbsM10000FileName.c_str());
    string incompatibleSitesFileName = vcfFileRoot + ".incompatibleSites";
    std::ofstream* incompatibleSitesFile = new std::ofstream(incompatibleSitesFileName.c_str());
    
    string line;
    // Read the location of the undetermined loci
    std::map<std::string, std::vector<std::string> > udlMap;
    std::map<std::string, int > scLengthMap;
    string thisScaffold;
    while (getline(*unlFile, line)) {
        if (line[0] == '>') {
            std::vector<std::string> scLvector = split(line, '\t');
            thisScaffold = scLvector[0].substr(1,string::npos);
            scLengthMap[thisScaffold] = atoi(scLvector[1].c_str());
        } else {
            udlMap[thisScaffold].push_back(line);
        }
    }
    
    // Now go through the VCF file
    std::vector<string> sampleNames;
    std::vector<string> fields;
    string currentScaffoldNum;
    int totalVariantNumber = 0;
    int numCombinations;
    std::map<string, std::map<string, std::vector<int> > > cbsTractLengths;
    std::map<string, std::map<string, std::vector<string> > > incompatibleSitesMap;
    std::vector<std::vector<string> > incompatibleSitesOneMatrix;
    std::map<string, int> currentCbsLengths;
    std::vector<std::vector<string> > toPrint;
    std::vector<std::vector<string> > toPrintMin10000;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            
        } else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
            // std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
            
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            assert(numSamples == sampleNames.size());
            
            numCombinations = choose((int)numSamples, 2);
            toPrint.resize(numCombinations); toPrintMin10000.resize(numCombinations);
            std::cerr << "Number of combinations: " << numCombinations << std::endl;
        } else {
            totalVariantNumber++;
            
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            int previousSNPpos;
            
            if (fields[0] != currentScaffoldNum) {
                if (currentScaffoldNum != "") {
                    // Summarise results for a scaffold
                    if (scLengthMap[currentScaffoldNum] >= opt::minScLength) {
                        std::vector<string> thisPair;
                        std::vector<string> thisPairMin10000;
                        int k = 0;
                        for (int i = 0; i != sampleNames.size(); i++) {
                            for (int j = i+1; j != sampleNames.size(); j++) {
                                std::vector<string> thisPair(numToStringVector(cbsTractLengths[sampleNames[i]+"+"+sampleNames[j]][currentScaffoldNum]));
                                std::vector<string> thisPairMin10000(numToStringMinFilterVector(cbsTractLengths[sampleNames[i]+"+"+sampleNames[j]][currentScaffoldNum], 10000));
                                if (currentScaffoldNum == "scaffold_0") {
                                    thisPair.insert(thisPair.begin(), sampleNames[i]+"+"+sampleNames[j]);
                                    thisPairMin10000.insert(thisPairMin10000.begin(), sampleNames[i]+"+"+sampleNames[j]);
                                }
                                toPrint[k].reserve(toPrint[k].size() + thisPair.size());
                                toPrint[k].insert(toPrint[k].end(), thisPair.begin(), thisPair.end());
                                toPrintMin10000[k].reserve(toPrintMin10000[k].size() + thisPairMin10000.size());
                                toPrintMin10000[k].insert(toPrintMin10000[k].end(), thisPairMin10000.begin(), thisPairMin10000.end());
                                k++;
                                //*cbsFile << sampleNames[i]+"+"+sampleNames[j] << "\t";
                                //print_vector(cbsTractLengths[sampleNames[i]+"+"+sampleNames[j]][currentScaffoldNum], *cbsFile);
                            }
                        }
                    } else {
                        std::cerr << "Sorting cbs tract length output... " << std::endl;
                        std::sort(toPrint.begin(), toPrint.end(), compareVectorsByLength);
                        std::sort(toPrintMin10000.begin(), toPrintMin10000.end(), compareVectorsByLength);
                        print_matrix(toPrint, *cbsFile);
                        print_matrix(toPrintMin10000, *cbsFileM10000);
                        int i = 0;
                        for (std::map<string, std::map<string, std::vector<string> > >::iterator it1 = incompatibleSitesMap.begin(); it1 != incompatibleSitesMap.end(); it1++) {
                            // *incompatibleSitesFile << it1->first;
                            for (std::map<string, std::vector<string> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                                incompatibleSitesOneMatrix[i].push_back(it1->first);
                                incompatibleSitesOneMatrix[i].push_back(it2->first);
                                incompatibleSitesOneMatrix[i].insert(incompatibleSitesOneMatrix[i].end(), it2->second.begin(), it2->second.end());
                                i++;
                                // *incompatibleSitesFile << it2->first;
                                // print_vector(it2->second, *incompatibleSitesFile);
                            }
                        }
                        std::cerr << "Sorting incompatible sites output... " << std::endl;
                        std::sort(incompatibleSitesOneMatrix.begin(), incompatibleSitesOneMatrix.end(), compareVectorsByLength);
                        print_matrix(incompatibleSitesOneMatrix, *incompatibleSitesFile);
                        exit(EXIT_SUCCESS);
                    }
                }
                currentScaffoldNum = fields[0];
                previousSNPpos = 1;
                for (int i = 0; i != sampleNames.size(); i++) {
                    for (int j = i+1; j != sampleNames.size(); j++) {
                        currentCbsLengths[sampleNames[i]+"+"+sampleNames[j]] = 0;
                    }
                }
                incompatibleSitesOneMatrix.resize(incompatibleSitesOneMatrix.size()+numCombinations);
                std::cerr << "Processing " << currentScaffoldNum << std::endl;
            }
            
            // Process variants
            if (info[0] != "INDEL") {  // Without indels
                Counts counts = getThisVariantCounts(fields);
                int thisSNPpos = atoi(fields[1].c_str());
                // Check compatibility for all pair of samples; 0 (hom_ref) is not compatible with 2 (hom_alt)
                for (int i = 0; i != sampleNames.size(); i++) {
                    for (int j = i+1; j != sampleNames.size(); j++) {
                        string thisCombination = sampleNames[i]+"+"+sampleNames[j];
                        if ((counts.individualsWithVariant[i] == 0 && counts.individualsWithVariant[j] == 2) ||
                            (counts.individualsWithVariant[i] == 2 && counts.individualsWithVariant[j] == 0)) {
                            cbsTractLengths[thisCombination][currentScaffoldNum].push_back(currentCbsLengths[thisCombination]);
                            incompatibleSitesMap[thisCombination][currentScaffoldNum].push_back(fields[1]);
                            currentCbsLengths[thisCombination] = 0;
                        } else {
                            currentCbsLengths[thisCombination] += (thisSNPpos - previousSNPpos);
                        }
                    }
                }
                previousSNPpos = thisSNPpos;
            }
            if (totalVariantNumber % 10000 == 0)
                std::cerr << totalVariantNumber << " variants processed..." << std::endl;
        }
    }
   
}


int cbsMain(int argc, char** argv) {
    parseCbsOptions(argc, argv);
    
    if (!opt::vcfFile.empty()) {
        std::cerr << "Calculating the lenghts of cbs tracts from: " << opt::vcfFile << std::endl;
        std::cerr << "and using information about scaffold lengths and the location of undetermined regions from: " << opt::undeterminedLociFile << std::endl;
        calculateCbs();
    } else if (!opt::genomeFile.empty()) {
        std::cerr << "Looking for undetermined regions (runs of Ns) in the reference genome: " << std::endl;
        std::cerr << opt::genomeFile << std::endl;
        findUndeterminedRegions();
    } else {
        assert(false);
    }
    
    return 0;
}

void parseCbsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 's': arg >> opt::sampleNameFile; break;
            case 'm': arg >> opt::minScLength; break;
            case OPT_CBS: arg >> opt::undeterminedLociFile; break;
            case OPT_PREPARE_GENOME: opt::bPrepareGenome = true; break;
            case 'h':
                std::cout << CBS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << CBS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    if (opt::bPrepareGenome)
        opt::genomeFile = argv[optind++];
    else if (!opt::undeterminedLociFile.empty()) {
        opt::vcfFile = argv[optind++];
    } else {
        std::cerr << "either the --cbs or the --prepare-genome option has to be selected" << std::endl;
        std::cout << CBS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
}
