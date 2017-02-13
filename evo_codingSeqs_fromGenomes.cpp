//
//  evo_codingSeqs_fromGenomes.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 06/02/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#include "evo_codingSeqs_fromGenomes.h"
#include "evo_protein_SegregatingSites.h"

#define SUBPROGRAM "SeqFromGenomes"

#define DEBUG 1

static const char *SEQFROMGENOME_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] <-g genome.fa | -l list_of_genomes_in_same_coordinates.txt> ANNOTATION.gffExtract\n"
"Get subsequences (usually genes) from whole genome sequences\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -g,   --genome FILE.fa                  a multiple alignment file (either -a or -l required)\n"
"       -l,   --listOfFiles LIST.txt            a list with multiple alignment filenames, one per line (either -a or -l required)\n"
"       --hapLabels=FILE_WITH_LABELS.txt        one per line, with as many records as there are haplotypes\n"
"       --outFolder=FOLDER                      put the output in FOLDER (default .)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hg:l:";

enum { OPT_NEW_LABELS, OPT_FOLDER };

static const struct option longopts[] = {
    { "genome",   required_argument, NULL, 'g' },
    { "listOfFiles",   required_argument, NULL, 'l' },
    { "newLabels", required_argument, NULL, OPT_NEW_LABELS },
    { "outFolder", required_argument, NULL, OPT_FOLDER },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string genomeFile = "";
    static string genomeListFile = "";
    static string newLabelFile = "";
    static string outFolder = "";
    static string geneFile = "";
}


int SeqFromGenomes(int argc, char** argv) {
    parseSeqFromGenomesOptions(argc, argv);
    // Load up the annotation file
    std::ifstream* geneFile = new std::ifstream(opt::geneFile.c_str());
    Annotation wgAnnotation(geneFile, false);
    std::vector<std::vector<string> > annotation;
    
    std::cerr << "Getting subsequences defined in " << opt::geneFile << std::endl;
    
    std::vector<string> allGenomeFiles; string AAFileName;
    std::vector<string> allNewLabels;
    if (opt::genomeFile != "") {
        allGenomeFiles.push_back(opt::genomeFile);
        std::cerr << "from the genome in: " << stripExtension(opt::genomeFile) << std::endl;
    } else {
        std::ifstream* alignmentList = new std::ifstream(opt::genomeListFile.c_str());
        std::cerr << "from the genomes listed in: " << opt::genomeListFile << std::endl;
        std::string line;
        while (getline(*alignmentList, line)) {
            allGenomeFiles.push_back(line);
        }
    }
    if (opt::newLabelFile != "") {
        std::cerr << "using labels from: " << opt::newLabelFile << std::endl;
        std::ifstream* newLabels = new std::ifstream(opt::newLabelFile.c_str());
        std::string line;
        while (getline(*newLabels, line)) {
            allNewLabels.push_back(line);
        }
    }
    
    // Loop over the genome files and load all into memory (so can handle only a few genomes):
    std::vector<std::map<std::string,string>> allGenomes;
    for (std::vector<std::string>::size_type i = 0; i != allGenomeFiles.size(); i++) {
        std::ifstream* genome = new std::ifstream(allGenomeFiles[i].c_str());
        std::map<std::string,string> allSeqs;
        // Load the whole genome
        std::string scaffoldSeq; scaffoldSeq.reserve(20000000); scaffoldSeq = "";
        std::string line; std::string thisScaffold;
        while (getline(*genome, line)) {
            if (line[0] == '>') {
                if (scaffoldSeq != "") { allSeqs[thisScaffold] = scaffoldSeq;}
                thisScaffold = line.substr(1); scaffoldSeq = "";
                continue;
            }
            scaffoldSeq.append(line);
        } genome->close();
        allSeqs[thisScaffold] = scaffoldSeq; // The final scaffold
        std::cerr << "Loaded sequences from: " << allGenomeFiles[i] << std::endl;
        allGenomes.push_back(allSeqs);
    }
    
    for(std::map<std::string,string>::iterator itGenome = allGenomes[0].begin(); itGenome != allGenomes[0].end(); ++itGenome) {
        //std::cerr << itGenome->first << std::endl;
        
        annotation = wgAnnotation.annotationMap[itGenome->first]; // Get annotation for this scaffold
        if (annotation.size() == 0) continue; // if we don't have any, just move to the next scaffold
        
        // check that this scaffold is present in all genomes:
        bool scaffoldInAllGenomes = true;
        for (std::vector<std::map<std::string,string>>::size_type i = 1; i != allGenomes.size(); i++) {
            if (allGenomes[i].count(itGenome->first) == 0)
                scaffoldInAllGenomes = false;
        }
        if (!scaffoldInAllGenomes) continue; // if not, just move on to the next scaffold
        
        //std::cerr << "Going through the annotation..." << std::endl;
        for (std::vector<std::vector<string> >::size_type k = 0; k != annotation.size(); k++) {
            std::vector<string> annotLineVec = split(annotation[k][0], '\t');
            string thisGeneName = annotLineVec[4];
            // std::cerr << "Gene:" << thisGeneName << std::endl;
            std::ofstream* geneOutFiles = new std::ofstream(thisGeneName.c_str());
            *geneOutFiles << ">" << stripExtension(allGenomeFiles[0]) << std::endl;
            string geneSequence = getReferenceForThisRegion(annotation[k], annotLineVec[3], itGenome->second);
            *geneOutFiles << geneSequence << std::endl;
            for (std::vector<std::map<std::string,string>>::size_type i = 1; i != allGenomes.size(); i++) {
                *geneOutFiles << ">" << stripExtension(allGenomeFiles[i]) << std::endl;
                string geneSequence = getReferenceForThisRegion(annotation[k], annotLineVec[3], allGenomes[i][itGenome->first]);
                *geneOutFiles << geneSequence << std::endl;
            } geneOutFiles->close();
            //std::cerr << "Got all sequences for: " << thisGeneName << std::endl;
        }
    }
    return 0;
}


void parseSeqFromGenomesOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'g': arg >> opt::genomeFile; break;
            case 'l': arg >> opt::genomeListFile; break;
            case OPT_NEW_LABELS: arg >> opt::newLabelFile; break;
            case OPT_FOLDER: arg >> opt::outFolder; break;
            case 'h':
                std::cout << SEQFROMGENOME_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::genomeFile == "" && opt::genomeListFile == "") {
        std::cerr << "Either -g or -l options must be specified\n";
        die = true;
    }
    
    if (opt::genomeFile != "" && opt::genomeListFile != "") {
        std::cerr << "The -g and -l options can't both be specified\n";
        die = true;
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
        std::cout << "\n" << SEQFROMGENOME_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    opt::geneFile = argv[optind++];
}
