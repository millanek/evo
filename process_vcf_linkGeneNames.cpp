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
"       --separateByCopyNumber=PREFIX           create separate files for gene IDs with 1-1, 1-N, N-1, and N-N relationships\n"
"                                               between the cichlid of interest and Zebrafish - the files will have the given PREFIX\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "ho:s:";
enum { OPT_V1, OPT_V2, OPT_NtoN, OPT_SEP_BY_COPY };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "v1",   required_argument, NULL, OPT_V1 },
    { "v2",   required_argument, NULL, OPT_V2 },
    { "separateByCopyNumber",   required_argument, NULL, OPT_SEP_BY_COPY },
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
    static string sepByCopyNumberPrefix= "";
}

inline int getSpeciesColumn(const string& species) {
    if (species == "mz") { return 0; }
    if (species == "pn") { return 1; }
    if (species == "ab") { return 2; }
    if (species == "nb") { return 3; }
    if (species == "on") { return 4; }
    return -1;
}


inline void attemptMappingUpdate(std::map<string,string>& cichlidHomologMap, const std::string& cichlidGene, const std::string& homologGene) {
    if (cichlidHomologMap.count(cichlidGene) == 0) {
        cichlidHomologMap[cichlidGene] = homologGene;
    } else {
        // Prefer zebrafish homologs as they have the best GO annotation
        if (cichlidHomologMap[cichlidGene].substr(0,6) == "ENSDAR") { return; }
        if (homologGene.substr(0,6) == "ENSDAR") { cichlidHomologMap[cichlidGene] = homologGene; }
        // Next best option is Medaka
        if (cichlidHomologMap[cichlidGene].substr(0,6) == "ENSORL") { return; }
        if (homologGene.substr(0,6) == "ENSORL") { cichlidHomologMap[cichlidGene] = homologGene; }
        // Next is stickleback:
        if (cichlidHomologMap[cichlidGene].substr(0,6) == "ENSGAC") { return; }
        if (homologGene.substr(0,6) == "ENSGAC") { cichlidHomologMap[cichlidGene] = homologGene; }
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
    
    // Load David Brawand's assignment of orthologs
    // Mapping from cichlid IDs to a zebrafish homolog (or medaka
    // stickleback, tetraodon, if zebrafish not available)
    std::map<string,string> cichlidHomolog;
    
    std::map<string,string> cichlidDanRerCopyNum;
    if (opt::v2orthologsFile != "") {
        std::cerr << "Reading the v2 full orthologs file:" << std::endl;
        std::ifstream* ocFile = new std::ifstream(opt::v2orthologsFile);
        while (getline(*ocFile, line)) {
            std::vector<string> orthVec = split(line, '\t');
            int c = getSpeciesColumn(opt::species);
            if (orthVec[c] != "NA") {
                if (orthVec[8] != "NA") { cichlidHomolog[orthVec[c]] = orthVec[8]; cichlidDanRerCopyNum[orthVec[c]] = "1-1"; } // Zebrafish
                else if (orthVec[5] != "NA") { cichlidHomolog[orthVec[c]] = orthVec[5]; } // Medaka
                else if (orthVec[7] != "NA") { cichlidHomolog[orthVec[c]] = orthVec[7]; } // Stickleback
                else if (orthVec[6] != "NA") { cichlidHomolog[orthVec[c]] = orthVec[6]; } // Tetraodon
                else { cichlidHomolog[orthVec[c]] = "novelCichlidGene"; }
            } else { continue; }
        } ocFile->close();
    }
    
    
    if (opt::v1orthologousClustersFile != "") {
        int copiesInCichlid = 0; int copiesInDanRer = 0;
        string cichlidGene = ""; string homologGene = "";
        
        std::cerr << "Reading the v1 orthologous cluster file: " << std::endl;
        std::ifstream* ocFile = new std::ifstream(opt::v1orthologousClustersFile);
        while (getline(*ocFile, line)) {
            std::vector<string> idAndNum = split(line, '\t');
            string thisLineGeneID = idAndNum[0];  int thisLineGeneClusterNumber = atoi(idAndNum[1].c_str());
            // Another line for the same cluster
            if (thisLineGeneClusterNumber == geneNum) {
                if (thisLineGeneID.substr(0,2) == opt::species) {
                    if (cichlidGene == "") {            // First copy in the cichlid species (e.g. mz)
                        cichlidGene = thisLineGeneID;
                    } else {                            // There is more than one copy in the cichlid
                        if (homologGene != "") {
                            if (copiesInDanRer <= 1 || opt::NtoN) {
                                attemptMappingUpdate(cichlidHomolog, cichlidGene, homologGene + "/" + numToString(copiesInCichlid));
                                if (copiesInDanRer == 1)
                                    cichlidDanRerCopyNum[cichlidGene] = "N-1";
                                else if (copiesInDanRer > 1)
                                    cichlidDanRerCopyNum[cichlidGene] = "N-N";
                            }
                            cichlidGene = thisLineGeneID;
                        }
                    }
                    copiesInCichlid++;
                } else if (thisLineGeneID.substr(0,6) == "ENSDAR") {
                    copiesInDanRer++;
                    if (homologGene == "") { homologGene = thisLineGeneID; }
                    else {
                        if (rand() < 0.5) homologGene = thisLineGeneID;  // 50% chance of using this zfish copy (hacky!!!)
                    }
                } else if (thisLineGeneID.substr(0,6) == "ENSGAC") {
                    if (homologGene == "") { homologGene = thisLineGeneID; }
                } else if (thisLineGeneID.substr(0,6) == "ENSORL") {
                    if (homologGene == "" || homologGene.substr(0,6) == "ENSGAC") { homologGene = thisLineGeneID; }
                } else if (thisLineGeneID.substr(0,6) == "ENSTNI") {
                    if (homologGene == "") { homologGene = thisLineGeneID; }
                }
                // std::cerr << atoi(idAndNum[1].c_str()) << "\t" << geneNum << std::endl;
            } else { // First line for a new cluster read
                // so first add the mapping for the previous cluster
                if (cichlidGene != "" && homologGene != "") {
                    assert(copiesInCichlid > 0);
                    if (copiesInDanRer == 1) {
                        if (copiesInCichlid == 1) {
                            cichlidDanRerCopyNum[cichlidGene] = "1-1";
                            attemptMappingUpdate(cichlidHomolog, cichlidGene,homologGene);
                        } else if (copiesInCichlid > 1) {
                            cichlidDanRerCopyNum[cichlidGene] = "N-1";
                            attemptMappingUpdate(cichlidHomolog, cichlidGene,homologGene + "/" + numToString(copiesInCichlid));
                        }
                    } else if (copiesInDanRer > 1) {
                        if (copiesInCichlid == 1) {
                            cichlidDanRerCopyNum[cichlidGene] = "1-N";
                            if (opt::NtoN)
                                attemptMappingUpdate(cichlidHomolog, cichlidGene,homologGene);
                        } else if (copiesInCichlid > 1) {
                            cichlidDanRerCopyNum[cichlidGene] = "N-N";
                            if (opt::NtoN)
                                attemptMappingUpdate(cichlidHomolog, cichlidGene,homologGene + "/" + numToString(copiesInCichlid));
                        }
                    } else {
                        if (copiesInCichlid == 1) {
                            attemptMappingUpdate(cichlidHomolog, cichlidGene,homologGene);
                        } else if (copiesInCichlid > 1) {
                            attemptMappingUpdate(cichlidHomolog, cichlidGene,homologGene + "/" + numToString(copiesInCichlid));
                        }
                    }
                }
                
                // then start looking through the next cluster
                cichlidGene = ""; homologGene = ""; copiesInDanRer = 0; copiesInCichlid = 0;
                geneNum = thisLineGeneClusterNumber;
                if (thisLineGeneID.substr(0,2) == opt::species) {
                    cichlidGene = thisLineGeneID;
                } else if (thisLineGeneID.substr(0,6) == "ENSDAR") {
                    copiesInDanRer++; homologGene = thisLineGeneID;
                } else if (thisLineGeneID.substr(0,6) == "ENSGAC") {
                    homologGene = thisLineGeneID;
                } else if (thisLineGeneID.substr(0,6) == "ENSORL") {
                    homologGene = thisLineGeneID;
                } else if (thisLineGeneID.substr(0,6) == "ENSTNI") {
                    homologGene = thisLineGeneID;
                }
            }
        } ocFile->close();
    }
    
    if (opt::sepByCopyNumberPrefix != "") {
        std::ofstream* OneOneFile = new std::ofstream(opt::sepByCopyNumberPrefix + "_1-1.txt");
        std::ofstream* NOneFile = new std::ofstream(opt::sepByCopyNumberPrefix + "_N-1.txt");
        std::ofstream* OneNFile = new std::ofstream(opt::sepByCopyNumberPrefix + "_1-N.txt");
        std::ofstream* NNFile = new std::ofstream(opt::sepByCopyNumberPrefix + "_N-N.txt");
        
        for (std::map<string, string>::iterator it = cichlidDanRerCopyNum.begin(); it != cichlidDanRerCopyNum.end(); it++) {
            if (it->second == "1-1") {
                *OneOneFile << it->first << std::endl;
            } else if (it->second == "N-1") {
                *NOneFile << it->first << std::endl;
            } else if (it->second == "1-N") {
                *OneNFile << it->first << std::endl;
            } else if (it->second == "N-N") {
                *NNFile << it->first << std::endl;
            }
        }
    
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
                print_vector(ensGene, std::cerr);
            }
           // std::cout << ensGene[0] << "\t" << ensGene[2] << std::endl;
        }
    }
    
  
    
    // Go through the gene prediction file and generate the final outputs
    std::ifstream* gpFile = new std::ifstream(opt::gpFile);
    int countNovel = 1; int countUnknown = 1; int countNotInEnsembl = 1;
    while (getline(*gpFile, line)) {
        std::vector<string> gpVec = split(line, '\t');
        if ( cichlidHomolog.count(gpVec[0]) == 1) {
            std::vector<string> ensembl = split(cichlidHomolog[gpVec[0]], '/');
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
            } else {
                std::cout << nameWdots << "\t" << "noOrthologAssigned" << "\t" << "0" << "\t" << opt::species + ".orthologNotInEnsembl." + numToString(countNotInEnsembl) << std::endl;
                *refLinkFile << opt::species + ".orthologNotInEnsembl." + numToString(countUnknown) << "\t" << "ortholog from Brawand data not foud in Ensembl v75" << "\t" << gpVec[0] << "\tNP_X\t77\t88\t" << "0" << "\t0" << std::endl;
                gpVec[11] = opt::species + ".orthologNotInEnsembl." + numToString(countNotInEnsembl);
                print_vector(gpVec, *gpOutFile);
                //std::cerr << ensembl[0] << std::endl;
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
            case OPT_SEP_BY_COPY: arg >> opt::sepByCopyNumberPrefix; break;
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
