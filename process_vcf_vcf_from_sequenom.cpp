//
//  process_vcf_vcf_from_sequenom.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 25/10/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_vcf_from_sequenom.h"
#define SUBPROGRAM "VCFfromSequenom"

static const char *VCF_FROM_SEQUENOM_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] SEQUENOM_OUTPUT.txt SEQUENOM_DESIGN.txt\n"
"Generate a VCF file from a sequenom output file\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1 };

static const char* shortopts = "hs:";

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "samples",   required_argument, NULL, 's' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string sequenomFile;
    static string sequenomDesignFile;
    static string runName;
}

int VCFfromSequenomMain(int argc, char** argv) {
    parseVCFfromSequenomOptions(argc, argv);
    string fileRoot = stripExtension(opt::sequenomFile);
    
    std::cerr << "Using genotype calls in: " << opt::sequenomFile << std::endl;
    std::cerr << "to generate a VCF file."<< std::endl;
    
    // Open connection to read from the vcf file
    std::ifstream* sequenomFile = new std::ifstream(opt::sequenomFile.c_str());
    std::ifstream* sequenomDesignFile = new std::ifstream(opt::sequenomDesignFile.c_str());
    string vcfFileName = fileRoot + opt::runName + ".vcf";
    std::ofstream* pVCF = new std::ofstream(vcfFileName.c_str());
    std::ofstream* pSampleNames = new std::ofstream("sample_names.txt");
    
    
    string line;
    std::map<string, std::vector<string> > varRefAlt;
    while (getline(*sequenomDesignFile, line)) {
        std::vector<std::string> fields = split(line, '\t');
        string variant = fields[0];
        fields = split(fields[1], '[');
        string baseBeforeVar = fields[0].substr(fields[0].length()-1,1);
        fields = split(fields[1], ']');
        fields = split(fields[0], '/');
        string ref = fields[0];
        string alt = fields[1];
        varRefAlt[variant].push_back(ref);
        varRefAlt[variant].push_back(alt);
        varRefAlt[variant].push_back(baseBeforeVar);
    }
    
    
    
    getline(*sequenomFile, line);
    std::string firstLocus;
    std::vector<std::string> samples;
    std::map<string, std::vector<string> > varGenotypes;
    std::vector<string> checkDuplicates;
    while (getline(*sequenomFile, line)) {
        std::vector<std::string> fields = split(line, '\t');
        string var = fields[5]; string gnt = fields[3]; string sample = fields[1];
        if (std::find(checkDuplicates.begin(),checkDuplicates.end(), sample+"\t"+var) != checkDuplicates.end()) {
            continue;
        } else {
            checkDuplicates.push_back(sample+"\t"+var);
        }
        //if (varGenotypes.find(fields[5]) != varGenotypes.end()) {
        string ref = varRefAlt[var][0]; string alt = varRefAlt[var][1]; string bbv = varRefAlt[var][2];
        if (varGenotypes[var].empty()) {
            std::vector<std::string> scafLoc = split(var, '_');
            if (ref == "-")
                varGenotypes[var].push_back(scafLoc[0]+"_"+scafLoc[1]+"\t"+scafLoc[2]+"\t"+"."+"\t"+bbv+"\t"+bbv+alt+"\t.\t.\t.\tGT");
            else if (alt == "-")
                varGenotypes[var].push_back(scafLoc[0]+"_"+scafLoc[1]+"\t"+scafLoc[2]+"\t"+"."+"\t"+bbv+ref+"\t"+bbv+"\t.\t.\t.\tGT");
            else
                varGenotypes[var].push_back(scafLoc[0]+"_"+scafLoc[1]+"\t"+scafLoc[2]+"\t"+"."+"\t"+ref+"\t"+alt+"\t.\t.\t.\tGT");
        }
        if ((gnt.length() == 2 && ref.length() == 1 && alt.length() == 1) || (gnt.find('.') != string::npos))
            varGenotypes[var].push_back("1/0");
        else if (fields[3] == varRefAlt[var][0])
            varGenotypes[var].push_back("0/0");
        else if (fields[3] == varRefAlt[var][1])
            varGenotypes[var].push_back("1/1");
        else if (gnt == "#NAME?" || gnt == "N")
            varGenotypes[var].push_back("./.");
        else
            assert(false); 
        if (firstLocus.empty()) {
            firstLocus = var;
            samples.push_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            samples.push_back(sample);
        } else if (firstLocus == fields[5]) {
            samples.push_back(sample);
        }
        
        fields = split(var, '_');
        //}
        
    }
    
    *pVCF << VCF_HEADER;
    print_vector(samples, *pVCF);
    for (std::map<string, std::vector<string> >::iterator it = varGenotypes.begin(); it != varGenotypes.end(); it++) {
        print_vector(it->second, *pVCF);
    }
    
    print_vector(samples, *pSampleNames, '\n');
    
    return 0;
}

void parseVCFfromSequenomOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'n': arg >> opt::runName; break;
            case '?': die = true; break;
            case 'h': std::cout << VCF_FROM_SEQUENOM_USAGE_MESSAGE;
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
        std::cout << "\n" << VCF_FROM_SEQUENOM_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::sequenomFile = argv[optind++];
    opt::sequenomDesignFile = argv[optind++];
}
