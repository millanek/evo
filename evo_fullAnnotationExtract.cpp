//
//  evo_fullAnnotationExtract.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 12/02/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#include "evo_fullAnnotationExtract.h"
#include <unordered_map>

#define SUBPROGRAM "AnnotationPreformat"

#define DEBUG 1

static const char *ANNOTEXTRACT_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] SINGLE_COVER_GENE_PRED_FILE.gp ANNOTATION.gtf\n"
"Extract CDS coordinates from ANNOTATION.gtf for transcripts defined in SINGLE_COVER_GENE_PRED_FILE.gp\n"
"SINGLE_COVER_GENE_PRED_FILE.gp can be obtained by running Jim Kent's 'genePredSingleCover'\n"
"the output of this program (name e.g. 'ANNOTATION.gtfExtract' can then be fed into 'evo getCodingSeq'\n"
"\n"
"       -h, --help                              display this help and exit\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "h";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string gpFile = "";
    static string annotationFile = "";
}


int AnnotationPreformat(int argc, char** argv) {
    parseAnnotExtractOptions(argc, argv);
    // Load up the annotation file
    std::ifstream* gpFile = new std::ifstream(opt::gpFile.c_str());
    std::unordered_map<string, bool> gpTranscripts;
    std::cerr << "Loading " << opt::gpFile << std::endl;
    std::string line;
    while (getline(*gpFile, line)) {
        gpTranscripts[split(line, '\t')[0]] = true;
    } gpFile->close();
    std::cerr << "Done" << std::endl;
    
    
    std::ifstream* annotFile = new std::ifstream(opt::annotationFile.c_str());
    std::cerr << "Going through " << opt::annotationFile << std::endl;
    while (getline(*annotFile, line)) {
        std::vector<string> annotVec = split(line, '\t');
        if (annotVec[2] == "CDS") { // only outputting CDS coordinates
            std::vector<string> descriptionVec = split(annotVec[8], ' ');
            string transcript = descriptionVec[3].substr(1).substr(0,descriptionVec[3].substr(1).length()-2);
            if (gpTranscripts.count(transcript) == 1) {
                string gene = descriptionVec[1].substr(1).substr(0,descriptionVec[1].substr(1).length()-2);
                string scaffold = annotVec[0];
                string start = annotVec[3];
                string end = annotVec[4];
                string direction = annotVec[6];
                std::cout << scaffold << "\t" << start << "\t" << end << "\t" << direction << "\t" << gene << std::endl;
            }
            //
        }
        // if (gpTranscripts.count() == 1) {
        // }
    } annotFile->close();
    std::cerr << "Done" << std::endl;
    return 0;
}


void parseAnnotExtractOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'h':
                std::cout << ANNOTEXTRACT_USAGE_MESSAGE;
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
        std::cout << "\n" << ANNOTEXTRACT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    opt::gpFile = argv[optind++];
    opt::annotationFile = argv[optind++];
}
