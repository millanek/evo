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
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] SINGLE_COVER_GENE_PRED_FILE.gp ci_orthologous_clusters.txt/full_orthologs ENSEMBL_Entrez_gene_link.txt\n"
"Use homology information from David Brawand to link his gene names with known zebrafish gene names\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -o, --out=RUN_NAME                      RUN_NAME will be a part of the output files' names\n"
"       -s, --species=SPECIES                   the cichlid species (mz,pn,ab,nb, or on)\n"
"       --v2                                    using v2 annotation (e.g. BROADMZ2,full_orthologs)\n"
"       --NtoN                                  include genes with 1-to-N and N-to-N relationships (for Gene Ontology analysis)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "ho:s:";
enum { OPT_V2, OPT_NtoN };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "v2",   no_argument, NULL, OPT_V2 },
    { "NtoN",   no_argument, NULL, OPT_NtoN },
    { "out",   no_argument, NULL, 'o' },
    { "species",   no_argument, NULL, 's' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    
    static string gpFile;
    static string orthologousClustersFile;
    static string ensGeneFile;
    static string ensEntrezFile;
    static string out = "";
    static bool v2 = false;
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
    bool multi = false;
    int copiesInCichlid = 1;
    string thisCichlid = ""; string thisDanRer = "";
    std::map<string,string> cichlidDanRer;
    std::ifstream* ocFile = new std::ifstream(opt::orthologousClustersFile);
    
    if (opt::v2) {
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
        }
    } else {
        while (getline(*ocFile, line)) {
            std::vector<string> idAndNum = split(line, '\t');
            if (atoi(idAndNum[1].c_str()) == geneNum) {
                if (idAndNum[0].substr(0,2) == opt::species) {
                    if (thisCichlid == "") { thisCichlid = idAndNum[0]; }
                    else { if (thisDanRer != "" && !multi) { cichlidDanRer[thisCichlid] = thisDanRer + "/" + numToString(copiesInCichlid); thisCichlid = idAndNum[0]; copiesInCichlid++;} }
                }
                if (idAndNum[0].substr(0,6) == "ENSDAR") {
                    if (thisDanRer == "") { thisDanRer = idAndNum[0]; }
                    else { if (!opt::NtoN) multi = true; else if (rand() < 0.5) thisDanRer = idAndNum[0];}
                }
                // std::cerr << atoi(idAndNum[1].c_str()) << "\t" << geneNum << std::endl;
            } else {
                if (thisCichlid != "" && thisDanRer != "" && !multi) { cichlidDanRer[thisCichlid] = thisDanRer; }
                thisCichlid = ""; thisDanRer = ""; multi = false; copiesInCichlid = 2;
                geneNum = atoi(idAndNum[1].c_str());
                if (idAndNum[0].substr(0,2) == opt::species) { thisCichlid = idAndNum[0]; }
                if (idAndNum[0].substr(0,6) == "ENSDAR") { thisDanRer = idAndNum[0]; }
            }
        }
    }
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
    
  
    
    
    std::ifstream* gpFile = new std::ifstream(opt::gpFile);
    int countNovel = 1;
    while (getline(*gpFile, line)) {
        std::vector<string> gpVec = split(line, '\t');
        if ( cichlidDanRer.find(gpVec[0]) != cichlidDanRer.end() ) {
            std::vector<string> ensembl = split(cichlidDanRer[gpVec[0]], '/');
            std::vector<string> myNameVec = split(gpVec[0], '.');
            gpVec[0] = myNameVec[0] + "_" + myNameVec[1] + "_" + myNameVec[2] + "_" + myNameVec[3];
            
            if ( ensGeneMap.find(ensembl[0]) != ensGeneMap.end() ) {
                if (ensembl.size() == 1) {
                    std::cout << gpVec[0] << "\t" << ensembl[0] << "\t" << ensGeneMap[ensembl[0]] << std::endl;
                    gpVec[11] = ensGeneMap[ensembl[0]];
                    print_vector(gpVec, *gpOutFile);
                    *refLinkFile << ensGeneMap[ensembl[0]] << "\t" << ensGeneDescriptionMap[ensembl[0]] << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << ensEntrezMap[ensembl[0]] << "\t0" << std::endl;
                    *fullBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    if (ensEntrezMap[ensembl[0]] != "0") {
                        *goBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    }
                } else {
                    std::cout << gpVec[0] << "\t" << ensembl[0] << "\t" << ensGeneMap[ensembl[0]] << "/" << ensembl[1] << std::endl;
                    gpVec[11] = ensGeneMap[ensembl[0]]+"/"+ensembl[1];
                    print_vector(gpVec, *gpOutFile);
                    *refLinkFile << ensGeneMap[ensembl[0]] << "/" << ensembl[1] << "\t" << ensGeneDescriptionMap[ensembl[0]] << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << ensEntrezMap[ensembl[0]] << "\t0" << std::endl;
                    *fullBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    if (ensEntrezMap[ensembl[0]] != "0") {
                        *goBedFile << gpVec[1] << "\t" << gpVec[3] << "\t" << gpVec[4] << "\t" << ensEntrezMap[ensembl[0]] << "\t0\t" << gpVec[2] << std::endl;
                    }
                }
            } else if (ensembl[0] == "novelCichlidGene") {
                std::cout << gpVec[0] << "\t" << ensembl[0] << "\t" << opt::species + ".novel." + numToString(countNovel) << std::endl;
                gpVec[11] = opt::species + ".novel." + numToString(countNovel);
                countNovel++;
                print_vector(gpVec, *gpOutFile);
                *refLinkFile << opt::species + ".novel." + numToString(countNovel) << "\t" << "novel gene found only in cichlids" << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << "0" << "\t0" << std::endl;
            }
            //std::cout << "hello" << std::endl;
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
            case OPT_V2: opt::v2 = true; break;
            case OPT_NtoN: opt::NtoN = true; break;
            case 'h':
                std::cout << LINKGN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
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
    opt::orthologousClustersFile = argv[optind++];
    if (argc - optind >= 1) {
        opt::ensGeneFile = argv[optind++];
    }
    
}