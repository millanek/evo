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
"       -r, --regulatory=BP_5',BP_3'            make a file with regulatory coordinates for each gene\n"
"                                               including BP_5' basepairs upstream (default 3000)\n"
"                                               and BP_3' basepairs downstream (default 1000)\n"
"                                               and introns\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "h";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "regulatory", optional_argument, NULL, 'r' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string gpFile = "";
    static string annotationFile = "";
    static bool regulatory = false;
    static int bp_5prime=3000;
    static int bp_3prime=1000;
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
    std::ofstream* regulatoryFile;
    std::ofstream* intronFile;
    std::ofstream* upstreamFile;
    std::ofstream* downstreamFile;
    if (opt::regulatory) {
        regulatoryFile = new std::ofstream(opt::annotationFile+"Extract_allRegulatory");
        intronFile = new std::ofstream(opt::annotationFile+"Extract_Intron");
        upstreamFile = new std::ofstream(opt::annotationFile+"Extract_Upstream");
        downstreamFile = new std::ofstream(opt::annotationFile+"Extract_Downstream");
    }
    string geneLastLine = ""; string scaffoldLastLine = ""; string startLastLine = ""; string endLastLine = "";
    string directionLastLine = "";
    std::cerr << "Going through " << opt::annotationFile << std::endl;
    while (getline(*annotFile, line)) {
        if (line[0] == '#') continue;
        std::vector<string> annotVec = split(line, '\t');
        if (annotVec[2] == "CDS") { // only outputting CDS coordinates
            std::vector<string> descriptionVec = split(annotVec[8], ' ');
            string transcript = "";
            if (descriptionVec[2] == "transcript_id") {
                transcript = descriptionVec[3].substr(1).substr(0,descriptionVec[3].substr(1).length()-2);
            } else if (descriptionVec[4] == "transcript_id") {
                transcript = descriptionVec[5].substr(1).substr(0,descriptionVec[5].substr(1).length()-2);
            }
            if (gpTranscripts.count(transcript) == 1) {
                string gene = descriptionVec[1].substr(1).substr(0,descriptionVec[1].substr(1).length()-2);
                string scaffold = annotVec[0];
                string start = annotVec[3];
                string end = annotVec[4];
                string direction = annotVec[6];
                std::cout << scaffold << "\t" << start << "\t" << end << "\t" << direction << "\t" << gene << std::endl;
                if (opt::regulatory) {
                    if (gene != geneLastLine) {
                        // First sort out the previous gene:
                        if (directionLastLine == "+") {
                            *downstreamFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_3prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "downstream" << std::endl;
                            *regulatoryFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_3prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "downstream" << std::endl;
                        }
                        if (directionLastLine == "-") {
                            *upstreamFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_5prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "upstream" << std::endl;
                            *regulatoryFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_5prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "upstream" << std::endl;
                        }
                        // Now sort out this gene:
                        if (direction == "+") {
                            int startUpstream = atoi(start.c_str())-opt::bp_5prime-1;
                            if (startUpstream < 0) startUpstream = 0;
                            int endUpstream = atoi(start.c_str())-1;
                            if (endUpstream > 0) {
                                *upstreamFile << scaffold << "\t" << startUpstream << "\t" << atoi(start.c_str())-1 << "\t" << direction << "\t" << gene << "\t" << "upstream" << std::endl;
                                *regulatoryFile << scaffold << "\t" << startUpstream << "\t" << atoi(start.c_str())-1 << "\t" << direction << "\t" << gene << "\t" << "upstream" << std::endl;
                            }
                        }
                        if (direction == "-") {
                            int startDownstream = atoi(start.c_str())-opt::bp_5prime-1;
                            if (startDownstream < 0) startDownstream = 0;
                            int endDownstream = atoi(start.c_str())-1;
                            if (endDownstream > 0) {
                                *downstreamFile << scaffold << "\t" << startDownstream << "\t" << atoi(start.c_str())-1 << "\t" << direction << "\t" << gene << "\t" << "downstream" << std::endl;
                                *regulatoryFile << scaffold << "\t" << startDownstream << "\t" << atoi(start.c_str())-1 << "\t" << direction << "\t" << gene << "\t" << "downstream" << std::endl;
                            }
                        }
                    }
                    if (gene == geneLastLine) {
                        int startIntron = atoi(endLastLine.c_str())+1;
                        int endIntron = atoi(start.c_str())-1;
                        if (endIntron > startIntron) {
                            *intronFile << scaffold << "\t" << startIntron << "\t" << endIntron << "\t" << direction << "\t" << gene << "\t" << "intron" << std::endl;
                            *regulatoryFile << scaffold << "\t" << startIntron << "\t" << endIntron << "\t" << direction << "\t" << gene << "\t" << "intron" << std::endl;
                        } else {
                            std::cerr << "WARNING: an intron for gene " << gene << " has negative length" << std::endl;
                        }
                    }
                }
                geneLastLine = gene; startLastLine = start; endLastLine = end; directionLastLine = direction;
                scaffoldLastLine = scaffold;
            }
            //
        }
        // if (gpTranscripts.count() == 1) {
        // }
    } annotFile->close();
    // and the final noncoding element:
    if (opt::regulatory) {
        if (directionLastLine == "+") {
            *downstreamFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_3prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "downstream" << std::endl;
            *regulatoryFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_3prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "downstream" << std::endl;
        }
        if (directionLastLine == "-") {
            *upstreamFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_5prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "upstream" << std::endl;
            *regulatoryFile << scaffoldLastLine << "\t" << atoi(endLastLine.c_str())+1 << "\t" << atoi(endLastLine.c_str())+opt::bp_5prime+1 << "\t" << directionLastLine << "\t" << geneLastLine << "\t" << "upstream" << std::endl;
        }
    }
    std::cerr << "Done" << std::endl;
    return 0;
}


void parseAnnotExtractOptions(int argc, char** argv) {
    bool die = false;
    std::string regulString = "";
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'r': opt::regulatory = true; arg >> regulString;
                break;
            case 'h':
                std::cout << ANNOTEXTRACT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (regulString != "") {
        std::vector<string> upDown = split(regulString, ',');
        if (upDown.size() != 2) {
            std::cerr << "Check the format of your --regulatory argument\n";
            die = true;
        } else {
            opt::bp_5prime = atoi(upDown[0].c_str()); opt::bp_3prime = atoi(upDown[1].c_str());
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
