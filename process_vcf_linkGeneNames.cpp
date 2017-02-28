//
//  process_vcf_linkGeneNames.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 25/03/2014.
//  Copyright (c) 2014 Milan Malinsky. All rights reserved.
//

#include "process_vcf_linkGeneNames.h"

#define SUBPROGRAM "linkGeneNames"

#define DEBUG 1

static const char *LINKGN_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] SINGLE_COVER_GENE_PRED_FILE.gp ENSEMBL_Entrez_gene_link.txt\n"
"Use homology information from David Brawand to link his gene names with known zebrafish gene names\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -o, --out=RUN_NAME                      RUN_NAME will be a part of the output files' names\n"
"       -s, --species=SPECIES                   the cichlid species (mz,pn,ab,nb, or on)\n"
"       --v1=ci_orthologous_clusters.txt        using ortholog assignment from v1 annotation\n"
"       --v2=full_orthologs                     using v2 annotation (e.g. BROADMZ2,full_orthologs)\n"
"                                               only --v1 can be used with BROADMZ1\n"
"                                               either --v1 or --v2 or both can be used with BROADMZ2\n"
"       --NtoN                                  include genes with 1-to-N and N-to-N relationships (for Gene Ontology analysis)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "ho:s:";
enum { OPT_V1, OPT_V2, OPT_NtoN };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "v1",   required_argument, NULL, OPT_V1 },
    { "v2",   required_argument, NULL, OPT_V2 },
    { "NtoN",   no_argument, NULL, OPT_NtoN },
    { "out",   no_argument, NULL, 'o' },
    { "species",   no_argument, NULL, 's' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    
    static string gpFile;
    static string v1orthologousClustersFile = "";
    static string ensGeneFile;
    static string ensEntrezFile;
    static string out = "";
    static string v2orthologsFile = "";
    static bool NtoN = false;
    static string species = "mz";
}

inline int getSpeciesColumn(const string& species) {
    if (species == "mz") { return 0; }
    if (species == "pn") { return 1; }
    if (species == "ab") { return 2; }
    if (species == "nb") { return 3; }
    if (species == "on") { return 4; }
    return -1;
}


inline void attemptMappingUpdate(std::map<string,string>& cichlidDanRerMap, const std::string& thisCichlid, const std::string& thisDanRer) {
    if (cichlidDanRerMap.count(thisCichlid) == 0) {
        cichlidDanRerMap[thisCichlid] = thisDanRer;
    } else {
        // Prefer zebrafish homologs as they have the best GO annotation
        if (cichlidDanRerMap[thisCichlid].substr(0,6) == "ENSDAR") { return; }
        if (thisDanRer.substr(0,6) == "ENSDAR") { cichlidDanRerMap[thisCichlid] = thisDanRer; }
        // Next best option is Medaka
        if (cichlidDanRerMap[thisCichlid].substr(0,6) == "ENSORL") { return; }
        if (thisDanRer.substr(0,6) == "ENSORL") { cichlidDanRerMap[thisCichlid] = thisDanRer; }
        // Next is stickleback:
        if (cichlidDanRerMap[thisCichlid].substr(0,6) == "ENSGAC") { return; }
        if (thisDanRer.substr(0,6) == "ENSGAC") { cichlidDanRerMap[thisCichlid] = thisDanRer; }
        // Tetraodon would be added only if there wasn't any other homolo before
    }
}

int linkGNMain(int argc, char** argv) {
    linkGNOptions(argc, argv);
    string gpFileRoot = stripExtension(opt::gpFile);
    std::ofstream* gpOutFile;
    std::ofstream* refLinkFile;
    std::ofstream* goBedFile;
    std::ofstream* fullBedFile;
    
    if (opt::NtoN) {
        goBedFile = new std::ofstream(gpFileRoot + opt::out + "_NtoN_GOBed.txt");
        fullBedFile = new std::ofstream(gpFileRoot + opt::out + "_NtoN_FullBed.txt");
        gpOutFile = new std::ofstream(gpFileRoot + opt::out + "_NtoN_RefGene.gp");
        refLinkFile = new std::ofstream(gpFileRoot + opt::out + "_NtoN_RefLink.gp");
    } else {
        goBedFile = new std::ofstream(gpFileRoot + opt::out + "_GOBed.txt");
        fullBedFile = new std::ofstream(gpFileRoot + opt::out + "_FullBed.txt");
        gpOutFile = new std::ofstream(gpFileRoot + opt::out + "_RefGene.gp");
        refLinkFile = new std::ofstream(gpFileRoot + opt::out + "_RefLink.gp");
    }
    
    string line;
    int geneNum = 1;
    bool multi = false; // If we don't want NtoN relationships, indicate that this has multiple copies in zfish
    int copiesInCichlid = 1;
    string thisCichlid = ""; string thisDanRer = "";
    
    // Load David Brawand's assignment of orthologs
    // Mapping from cichlid IDs to a zebrafish ortholog (or medaka
    // stickleback, tetraodon, if zebrafish not available)
    std::map<string,string> cichlidDanRer;
    if (opt::v2orthologsFile != "") {
        std::cerr << "Reading the v2 full orthologs file:" << std::endl;
        std::ifstream* ocFile = new std::ifstream(opt::v2orthologsFile);
        while (getline(*ocFile, line)) {
            std::vector<string> orthVec = split(line, '\t');
            int c = getSpeciesColumn(opt::species);
            if (orthVec[c] != "NA") {
                if (orthVec[8] != "NA") { cichlidDanRer[orthVec[c]] = orthVec[8]; } // Zebrafish
                else if (orthVec[5] != "NA") { cichlidDanRer[orthVec[c]] = orthVec[5]; } // Medaka
                else if (orthVec[7] != "NA") { cichlidDanRer[orthVec[c]] = orthVec[7]; } // Stickleback
                else if (orthVec[6] != "NA") { cichlidDanRer[orthVec[c]] = orthVec[6]; } // Tetraodon
                else { cichlidDanRer[orthVec[c]] = "novelCichlidGene"; }
            } else { continue; }
        } ocFile->close();
    }
    
    if (opt::v1orthologousClustersFile != "") {
        std::cerr << "Reading the v1 orthologous cluster file: " << std::endl;
        std::ifstream* ocFile = new std::ifstream(opt::v1orthologousClustersFile);
        while (getline(*ocFile, line)) {
            std::vector<string> idAndNum = split(line, '\t');
            string thisLineGeneID = idAndNum[0];  int thisLineGeneClusterNumber = atoi(idAndNum[1].c_str());
            // Another line for the same cluster
            if (thisLineGeneClusterNumber == geneNum) {
                if (thisLineGeneID.substr(0,2) == opt::species) {
                    // First copy in the cichlid species (e.g. mz)
                    if (thisCichlid == "") { thisCichlid = thisLineGeneID; }
                    else { // There is more than one copy in the cichlid
                        if (thisDanRer != "" && !multi) {
                            attemptMappingUpdate(cichlidDanRer, thisCichlid, thisDanRer + "/" + numToString(copiesInCichlid));
                            thisCichlid = thisLineGeneID; copiesInCichlid++;
                        }
                    }
                }
                if (thisLineGeneID.substr(0,6) == "ENSDAR") {
                    if (thisDanRer == "") { thisDanRer = thisLineGeneID; }
                    else {
                        if (!opt::NtoN) multi = true; // Indicate that this gene has multiple copies in zfish
                        else if (rand() < 0.5) thisDanRer = thisLineGeneID; // 50% chance of using this zfish copy
                    }
                } else if (thisLineGeneID.substr(0,6) == "ENSGAC") {
                    if (thisDanRer == "") { thisDanRer = thisLineGeneID; }
                } else if (thisLineGeneID.substr(0,6) == "ENSORL") {
                    if (thisDanRer == "" || thisDanRer.substr(0,6) == "ENSGAC") { thisDanRer = thisLineGeneID; }
                } else if (thisLineGeneID.substr(0,6) == "ENSTNI") {
                    if (thisDanRer == "") { thisDanRer = thisLineGeneID; }
                }
                // std::cerr << atoi(idAndNum[1].c_str()) << "\t" << geneNum << std::endl;
            } else { // First line for a new cluster read
                // so first add the mapping for the previous cluster
                if (thisCichlid != "" && thisDanRer != "" && !multi) {
                    attemptMappingUpdate(cichlidDanRer, thisCichlid,thisDanRer);
                }
                thisCichlid = ""; thisDanRer = ""; multi = false; copiesInCichlid = 2;
                geneNum = thisLineGeneClusterNumber;
                if (thisLineGeneID.substr(0,2) == opt::species) { thisCichlid = thisLineGeneID; }
                if (thisLineGeneID.substr(0,6) == "ENSDAR") { thisDanRer = thisLineGeneID; }
            }
        } ocFile->close();
    }
    
    // Load gene names and descriptions from ENSEMBL
    std::map<string,string> ensGeneMap;
    std::map<string,string> ensGeneDescriptionMap;
    std::map<string,string> ensEntrezMap;
    if (!opt::ensGeneFile.empty()) {
        std::ifstream* egFile = new std::ifstream(opt::ensGeneFile);
        while (getline(*egFile, line)) {
            std::vector<string> ensGene = split(line, '\t');
            if (ensGene.size() == 4) {
                ensGeneMap[ensGene[0]] = ensGene[3];
                ensGeneDescriptionMap[ensGene[0]] = ensGene[2];
                // Sometimes there are two Entrez records for one Ensembl gene, the first Entrez record tends to be the more informative one
                if ( ensEntrezMap.find(ensGene[0]) == ensEntrezMap.end() ) {
                    if (ensGene[1] != "") {ensEntrezMap[ensGene[0]] = ensGene[1]; }
                    else { ensEntrezMap[ensGene[0]] = "0"; }
                }
            } else if (ensGene.size() == 3) {
                ensGeneMap[ensGene[0]] = "NA";
                if (ensGene[2] != "") { ensGeneDescriptionMap[ensGene[0]] = ensGene[2]; }
                else { ensGeneDescriptionMap[ensGene[0]] = "no description: " + ensGene[0]; }
                // Sometimes there are two Entrez records for one Ensembl gene, the first Entrez record tends to be the more informative one
                if ( ensEntrezMap.find(ensGene[0]) == ensEntrezMap.end() ) {
                    if (ensGene[1] != "") {ensEntrezMap[ensGene[0]] = ensGene[1]; }
                    else { ensEntrezMap[ensGene[0]] = "0"; }
                }
            } else {
                //std::cerr << ensGene.size() << std::endl;
                print_vector_stream(ensGene, std::cerr);
            }
           // std::cout << ensGene[0] << "\t" << ensGene[2] << std::endl;
        }
    }
    
  
    
    // Go through the gene prediction file and generate the final outputs
    std::ifstream* gpFile = new std::ifstream(opt::gpFile);
    int countNovel = 1; int countUnknown = 1;
    while (getline(*gpFile, line)) {
        std::vector<string> gpVec = split(line, '\t');
        if ( cichlidDanRer.count(gpVec[0]) == 1) {
            std::vector<string> ensembl = split(cichlidDanRer[gpVec[0]], '/');
            std::vector<string> myNameVec = split(gpVec[0], '.');
            std::string nameWdots = gpVec[0];
            gpVec[0] = myNameVec[0] + "_" + myNameVec[1] + "_" + myNameVec[2] + "_" + myNameVec[3];
            
            if ( ensGeneMap.find(ensembl[0]) != ensGeneMap.end() ) {
                if (ensembl.size() == 1) {
                    std::cout << nameWdots << "\t" << ensembl[0] << "\t" << ensEntrezMap[ensembl[0]] << "\t" << ensGeneMap[ensembl[0]] << std::endl;
                    gpVec[11] = ensGeneMap[ensembl[0]];
                    print_vector(gpVec, *gpOutFile);
                    *refLinkFile << ensGeneMap[ensembl[0]] << "\t" << ensGeneDescriptionMap[ensembl[0]] << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << ensEntrezMap[ensembl[0]] << "\t0" << std::endl;
                    *fullBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    if (ensEntrezMap[ensembl[0]] != "0") {
                        *goBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    }
                } else {
                    std::cout << nameWdots << "\t" << ensembl[0] << "\t" << ensEntrezMap[ensembl[0]] << "\t" << ensGeneMap[ensembl[0]] << "/" << ensembl[1] << std::endl;
                    gpVec[11] = ensGeneMap[ensembl[0]]+"/"+ensembl[1];
                    print_vector(gpVec, *gpOutFile);
                    *refLinkFile << ensGeneMap[ensembl[0]] << "/" << ensembl[1] << "\t" << ensGeneDescriptionMap[ensembl[0]] << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << ensEntrezMap[ensembl[0]] << "\t0" << std::endl;
                    *fullBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    if (ensEntrezMap[ensembl[0]] != "0") {
                        *goBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    }
                }
            } else if (ensembl[0] == "novelCichlidGene") {
                std::cout << nameWdots << "\t" << ensembl[0] << "\t0" << "\t" << opt::species + ".novel." + numToString(countNovel) << std::endl;
                gpVec[11] = opt::species + ".novel." + numToString(countNovel);
                countNovel++;
                print_vector(gpVec, *gpOutFile);
                *refLinkFile << opt::species + ".novel." + numToString(countNovel) << "\t" << "novel gene found only in cichlids" << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << "0" << "\t0" << std::endl;
            }
            //std::cout << "hello" << std::endl;
        } else {
            std::vector<string> myNameVec = split(gpVec[0], '.');
            std::string nameWdots = gpVec[0] + ".1";
            gpVec[0] = myNameVec[0] + "_" + myNameVec[1] + "_" + myNameVec[2] + "_" + myNameVec[3];
            std::cout << nameWdots << "\t" << "noOrthologAssigned" << "\t" << "0" << "\t" << opt::species + ".unknown." + numToString(countUnknown) << std::endl;
            *refLinkFile << opt::species + ".unknown." + numToString(countUnknown) << "\t" << "unknown - no ortholog from Brawand data" << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << "0" << "\t0" << std::endl;
            gpVec[11] = opt::species + ".unknown." + numToString(countUnknown);
            print_vector(gpVec, *gpOutFile);
            countUnknown++;
        }
    }
    
    return 0;
}

void linkGNOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'o': arg >> opt::out; break;
            case 's': arg >> opt::species; break;
            case OPT_V1: arg >> opt::v1orthologousClustersFile; break;
            case OPT_V2: arg >> opt::v2orthologsFile; break;
            case OPT_NtoN: opt::NtoN = true; break;
            case 'h':
                std::cout << LINKGN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::v1orthologousClustersFile == "" && opt::v2orthologsFile == "") {
        std::cerr << "either --v1 or --v2 ortholog assigment must be supplied\n";
        die = true;
    }
    
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (opt::out != "") opt::out = "_" + opt::out;
    
    if (die) {
        std::cout << "\n" << LINKGN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::gpFile = argv[optind++];
    if (argc - optind >= 1) {
        opt::ensGeneFile = argv[optind++];
    }
    
}
