//
//  evo_diversityFromHaps.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 04/03/2022.
//  Copyright © 2022 Milan Malinsky. All rights reserved.
//

#include "evo_diversityFromHaps.h"


#define SUBPROGRAM "RegionsPiGeneral"

#define DEBUG 1

static const char *REGION_PI_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] regions.bed variants.vcf\n"
"Calculate Dxy (and possibly other statistics) for regions defined in the .bed file\n"
"\n"
"       -h, --help                          display this help and exit\n"
"       -o, --outfile=FILE                  (optional; default=) output file name\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "ho:";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "outfile", required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string bedFile = "";
    static string vcfFile = "";
    static string outFile = "";
    static bool bElements = false;
}


int regionPiMain(int argc, char** argv) { 
    parseRegionPiOptions(argc, argv);
    
    // Load up the intervals file
    std::ifstream* bedFile = new std::ifstream(opt::bedFile.c_str());
    std::map<int, string> linearToGenomeMap;
    
    std::cerr << "Loading coordinates from the bed file " << opt::bedFile << std::endl;
    LinkedCoordsBed coords = LinkedCoordsBed(bedFile);          // The bed file needs element names in the fourth column
    std::cerr << "Done" << std::endl;
    
    // Load up the VCF file:
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    string line; // std::unordered_map<string, MultiallelicCounts> vcfCountsMap;
    std::unordered_map<string, double> vcfPiMap;
    std::unordered_map<string, double> vcfHetMap;
    std::unordered_map<string, double> vcfInbreedingPvalMap;
    int totalVariantNumber = 0;
    
    std::cerr << "Loading variants from the VCF file " << opt::vcfFile << std::endl;
    std::vector<std::string> fields; int numSamples = -1;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
        
        } else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            numSamples = (int)sampleNames.size();
        } else {
            totalVariantNumber++;
            fields = split(line, '\t');
            string scaffold = fields[0]; string pos = fields[1]; string altAllele = fields[4];
            std::vector<string> altAlleles = split(altAllele, ','); int starPos = -1;
            std::vector<std::string>::iterator it = std::find(altAlleles.begin(), altAlleles.end(), "*");
            if (it != altAlleles.end()) { starPos = (int) std::distance(altAlleles.begin(), it); }
          /*  if (starPos != -1) {
                std::cerr << "starPos: " << starPos << std::endl;
                std::cerr << "altAllele: " << altAllele << std::endl;
            } */
            
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            
            MultiallelicCounts* counts = new MultiallelicCounts(numSamples,starPos);
            counts->getMultiallelicCounts(genotypes);
            if (totalVariantNumber % 10000 == 0) {
                std::cerr << "Variant N: " << totalVariantNumber << std::endl;
            }
            double thisVariantPi = counts->getPiThisVariant();
            double thisVariantHet = counts->getHeterozygosityThisVariant();
            //vcfCountsMap[scaffold+"\t"+pos] = *counts;
            vcfPiMap[scaffold+"\t"+pos] = thisVariantPi;
            vcfHetMap[scaffold+"\t"+pos] = thisVariantHet;
            delete counts;
        }
    }
    std::cerr << "Done" << std::endl;
    

    string outPiFileName;
    if (opt::outFile != "") {
        outPiFileName = opt::outFile;
    } else {
        size_t suffixPos = opt::bedFile.find_last_of('.');
        if (suffixPos != std::string::npos) {
            outPiFileName = stripExtension(opt::bedFile) + "_PiVals_perElement.txt";
        } else {
            outPiFileName = opt::bedFile + "_PiVals_perElement.txt";
        }
    }
    std::ofstream* outPiFile = new std::ofstream(outPiFileName.c_str());
    
    std::cerr << "Now calculating Pi and Heterozygosity values per element to be written out into: " << outPiFileName << std::endl;
    std::vector<double> elementPiValues = coords.getMeanPerElement(vcfPiMap);
    std::vector<double> elementHetValues = coords.getMeanPerElement(vcfHetMap);
    std::cerr << "Calculations done" << std::endl;
    
    std::vector<string> elementNames = coords.getElementNames();
    std::vector<std::vector<string> > elementOuterBounds = coords.getElementOuterBoundaries();
    for (int i = 0; i < (int)elementPiValues.size(); i++) {
        *outPiFile << elementOuterBounds[i][0] << "\t" << elementOuterBounds[i][1] << "\t" << elementOuterBounds[i][2] << "\t" << elementNames[i] << "\t" << elementPiValues[i] << "\t" << elementHetValues[i] << std::endl;
    }
    
    std::cerr << "Done" << std::endl;
    return 0;
}


void parseRegionPiOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::outFile; break;
            case 'h':
                std::cout << REGION_PI_USAGE_MESSAGE;
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
        std::cout << "\n" << REGION_PI_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    opt::bedFile = argv[optind++];
    opt::vcfFile = argv[optind++];
}
