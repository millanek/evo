//
//  evo_Dmin_combine.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 13/07/2018.
//  Copyright © 2018 Milan Malinsky. All rights reserved.
//

#include "evo_Dmin_combine.h"


#define SUBPROGRAM "DminCombine"

#define DEBUG 1

static const char *DMINCOMBINE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] DminFile1 DminFile2 DminFile3 ....\n"
"Combine the BBAA, ABBA, and BABA counts from multiple files (e.g per-chromosome) and output the overall Dmin stats\n"
"also the D stats for the trio arrangement where the BBAA is the most common pattern\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_AA_EQ_O };

static const char* shortopts = "hn:";

static const int JK_WINDOW = 20000;

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static std::vector<string> dminFiles;
    static string runName = "combined";
}


int DminCombineMain(int argc, char** argv) {
    parseDminCombineOptions(argc, argv);
    string line; // for reading the input files
    
    
    std::vector<std::istream*> dminstdErrFiles; std::vector<std::istream*> dminBBAAscoreFiles;
    for (int i = 0; i < opt::dminFiles.size(); i++) {
        std::istream* dminBBAAscoreFile = createReader((opt::dminFiles[i] + "_combine.txt").c_str());
        dminBBAAscoreFiles.push_back(dminBBAAscoreFile);
        std::istream* dminstdErrFile;
        if (file_exists(opt::dminFiles[i] + "_combine_stderr.txt")) {
            dminstdErrFile = createReader((opt::dminFiles[i] + "_combine_stderr.txt").c_str());
        } else if(file_exists(opt::dminFiles[i] + "_combine_stderr.txt.gz")) {
            dminstdErrFile = createReader((opt::dminFiles[i] + "_combine_stderr.txt.gz").c_str());
        } else {
            std::cerr << "Can't fine the file: " << opt::dminFiles[i] + "_combine_stderr.txt" << " or " << opt::dminFiles[i] + "_combine_stderr.txt.gz" << std::endl;
        }
        dminstdErrFiles.push_back(dminstdErrFile);
        std::cerr << "Reading file " << opt::dminFiles[i] << std::endl;
    }
    
    
    // Now get the standard error values
    std::ofstream* outFileBBAA = new std::ofstream(opt::runName + "_BBAA.txt"); std::ofstream* outFileDmin = new std::ofstream(opt::runName + "_Dmin.txt");
    std::vector<double> BBAA_local_Ds; std::vector<double> ABBA_local_Ds; std::vector<double> BABA_local_Ds;
    string s1; string s2; string s3;
    double BBAAtotal = 0; double ABBAtotal = 0; double BABAtotal = 0;
    bool allDone = false;
    do {
        for (int i = 0; i < dminBBAAscoreFiles.size(); i++) {
            if (getline(*dminBBAAscoreFiles[i], line)) {
                std::vector<string> patternCounts = split(line, '\t');
                assert(patternCounts.size() == 6);
                if (i == 0) {
                    s1 = patternCounts[0];
                    s2 = patternCounts[1];
                    s3 = patternCounts[2];
                } else {
                    assert(s1 == patternCounts[0]); assert(s2 == patternCounts[1]); assert(s3 == patternCounts[2]);
                }
                double BBAA = stringToDouble(patternCounts[3]);
                double BABA = stringToDouble(patternCounts[4]);
                double ABBA = stringToDouble(patternCounts[5]);
                BBAAtotal += BBAA; ABBAtotal += ABBA; BABAtotal += BABA;
            }
        }
        //std::cerr << "ABBAtotal = " << ABBAtotal << std::endl;
        //std::cerr << "BABAtotal = " << BABAtotal << std::endl;
        double Dnum1 = ABBAtotal - BABAtotal; // assert(Dnum1 == Dnums[i][0]);
        double Dnum2 = ABBAtotal - BBAAtotal; // assert(Dnum2 == Dnums[i][1]);
        double Dnum3 = BBAAtotal - BABAtotal; // assert(Dnum3 == Dnums[i][2]);
        double Ddenom1 = ABBAtotal + BABAtotal; // assert(Ddenom1 == Ddenoms[i][0]);
        double Ddenom2 = ABBAtotal + BBAAtotal; // assert(Ddenom2 == Ddenoms[i][1]);
        double Ddenom3 = BBAAtotal + BABAtotal; // assert(Ddenom3 == Ddenoms[i][2]);
        //std::cerr << "Dnum1 = " << Dnum1 << std::endl;
        //std::cerr << "Ddenom1 = " << Ddenom1 << std::endl;
        double D1 = Dnum1/Ddenom1; double D2 = Dnum2/Ddenom2; double D3 = Dnum3/Ddenom3;
        //std::cerr << "D1 = " << D1 << std::endl;
        for (int i = 0; i < dminstdErrFiles.size(); i++) {
            if (getline(*dminstdErrFiles[i], line)) {
                std::vector<string> localDs = split(line, '\t');
                assert(localDs.size() == 3 || localDs.size() == 0);
                if (localDs.size() == 3) {
                    std::vector<string> BBAA_D_strings = split(localDs[0], ',');
                    std::vector<string> BABA_D_strings = split(localDs[1], ',');
                    std::vector<string> ABBA_D_strings = split(localDs[2], ',');
                    for (int j = 0; j < BBAA_D_strings.size(); j++) {
                        //std::cerr << "BBAA_D_strings[j] = " << BBAA_D_strings[j] << std::endl;
                        double thisBBAA_localD = stringToDouble(BBAA_D_strings[j]);
                        if (!std::isnan(thisBBAA_localD)) BBAA_local_Ds.push_back(thisBBAA_localD);
                       // std::cerr << "BABA_D_strings[j] = " << BABA_D_strings[j] << std::endl;
                        double thisBABA_localD = stringToDouble(BABA_D_strings[j]);
                        if (!std::isnan(thisBABA_localD)) BABA_local_Ds.push_back(thisBABA_localD);
                        //std::cerr << "ABBA_D_strings[j] = " << ABBA_D_strings[j] << std::endl;
                        double thisABBA_localD = stringToDouble(ABBA_D_strings[j]);
                        if (!std::isnan(thisABBA_localD)) ABBA_local_Ds.push_back(thisABBA_localD);
                    }
                }
            } else { // all lines have been processed
                allDone = true; break;
            }
        }
       // std::cerr << "D1 = " << D1 << std::endl;
        //print_vector_stream(BBAA_local_Ds, std::cerr);
        double BBAAstdErr = jackknive_std_err(BBAA_local_Ds);
        //print_vector_stream(BABA_local_Ds, std::cerr);
        double BABAstdErr = jackknive_std_err(BABA_local_Ds);
        //print_vector_stream(ABBA_local_Ds, std::cerr);
        //std::cerr << "D1 = " << D1 << std::endl;
        double ABBAstdErr = jackknive_std_err(ABBA_local_Ds);
        //std::cerr << "D1 = " << D1 << std::endl;
        //std::cerr << "BBAAstdErr" << BBAAstdErr << std::endl;
        double D1_Z = fabs(D1)/BBAAstdErr; double D2_Z = fabs(D2)/BABAstdErr;
        double D3_Z = fabs(D3)/ABBAstdErr;
        //std::cerr << "D1_Z = " << D1_Z << std::endl;
        if (s1 == "Altcal" && s2 == "Altshe" && s3 == "Asplep") {
            std::cerr << "D1_Z = " << D1_Z << std::endl;
            std::cerr << "BBAAstdErr = " << BBAAstdErr << std::endl;
            print_vector_stream(BBAA_local_Ds, std::cerr);
            std::cerr << "ABBAstdErr = " << ABBAstdErr << std::endl;
            std::cerr << "BABAstdErr = " << BABAstdErr << std::endl;
            std::cerr << "ABBAtotal = " << ABBAtotal << std::endl;
            std::cerr << "BABAtotal = " << BABAtotal << std::endl;
            std::cerr << "BBAAtotal = " << BBAAtotal << std::endl;
        }

        
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotal >= BABAtotal && BBAAtotal >= ABBAtotal) {
            if (D1 >= 0)
                *outFileBBAA << s1 << "\t" << s2 << "\t" << s3;
            else
                *outFileBBAA << s2 << "\t" << s1 << "\t" << s3;
            *outFileBBAA << "\t" << fabs(D1) << "\t" << D1_Z << "\t";
            *outFileBBAA << BBAAtotal << "\t" << BABAtotal << "\t" << ABBAtotal << std::endl;
        } else if (BABAtotal >= BBAAtotal && BABAtotal >= ABBAtotal) {
            if (D2 >= 0)
                *outFileBBAA << s1 << "\t" << s3 << "\t" << s2;
            else
                *outFileBBAA << s3 << "\t" << s1 << "\t" << s2;
            *outFileBBAA << "\t" << fabs(D2) << "\t" << D2_Z << "\t";
            *outFileBBAA << BABAtotal << "\t" << BBAAtotal << "\t" << ABBAtotal << std::endl;
        } else if (ABBAtotal >= BBAAtotal && ABBAtotal >= BABAtotal) {
            if (D3 >= 0)
                *outFileBBAA << s3 << "\t" << s2 << "\t" << s1;
            else
                *outFileBBAA << s2 << "\t" << s3 << "\t" << s1;
            *outFileBBAA << "\t" << fabs(D3) << "\t" << D3_Z << "\t";
            *outFileBBAA << ABBAtotal << "\t" << BABAtotal << "\t" << BBAAtotal << std::endl;
        }
        
        // Find Dmin:
        if (fabs(D1) <= fabs(D2) && fabs(D1) <= fabs(D3)) { // (P3 == S3)
            if (D1 >= 0)
                *outFileDmin << s1 << "\t" << s2 << "\t" << s3 << "\t" << D1 << "\t" << D1_Z << "\t" << std::endl;
            else
                *outFileDmin << s1 << "\t" << s2 << "\t" << s3 << "\t" << fabs(D1) << "\t" << D1_Z << "\t"<< std::endl;
        } else if (fabs(D2) <= fabs(D1) && fabs(D2) <= fabs(D3)) { // (P3 == S2)
            if (D2 >= 0)
                *outFileDmin << s1 << "\t" << s3 << "\t" << s2 << "\t" << D2 << "\t" << D2_Z << "\t"<< std::endl;
            else
                *outFileDmin << s3 << "\t" << s1 << "\t" << s2 << "\t" << fabs(D2) << "\t" << D2_Z << "\t"<< std::endl;
        } else if (fabs(D3) <= fabs(D1) && fabs(D3) <= fabs(D2)) { // (P3 == S1)
            if (D3 >= 0)
                *outFileDmin << s3 << "\t" << s2 << "\t" << s1 << "\t" << D3 << "\t" << D3_Z << "\t"<< std::endl;
            else
                *outFileDmin << s2 << "\t" << s3 << "\t" << s1 << "\t" << fabs(D3) << "\t" << D3_Z << "\t" << std::endl;;
        }
        
        BBAA_local_Ds.clear(); ABBA_local_Ds.clear(); BABA_local_Ds.clear();
        BBAAtotal = 0; ABBAtotal = 0; BABAtotal = 0;
    } while(!allDone);
  
    
    
    /*
    for (int i = 0; i != trios.size(); i++) { //
        // Get the standard error values:
        double D1stdErr = jackknive_std_err(regionDs[i][0]); double D2stdErr = jackknive_std_err(regionDs[i][1]);
        double D3stdErr = jackknive_std_err(regionDs[i][2]);
        // Get the D values
        //Dnums[i][0] = ABBAtotals[i] - BABAtotals[i]; Dnums[i][1] = ABBAtotals[i] - BBAAtotals[i]; Dnums[i][2] = BBAAtotals[i] - BABAtotals[i];
        double Dnum1 = ABBAtotals[i] - BABAtotals[i]; // assert(Dnum1 == Dnums[i][0]);
        double Dnum2 = ABBAtotals[i] - BBAAtotals[i]; // assert(Dnum2 == Dnums[i][1]);
        double Dnum3 = BBAAtotals[i] - BABAtotals[i]; // assert(Dnum3 == Dnums[i][2]);
        // Ddenoms[i][0] = ABBAtotals[i] + BABAtotals[i]; Ddenoms[i][1] = ABBAtotals[i] + BBAAtotals[i]; Ddenoms[i][2] = BBAAtotals[i] + BABAtotals[i];
        double Ddenom1 = ABBAtotals[i] + BABAtotals[i]; // assert(Ddenom1 == Ddenoms[i][0]);
        double Ddenom2 = ABBAtotals[i] + BBAAtotals[i]; // assert(Ddenom2 == Ddenoms[i][1]);
        double Ddenom3 = BBAAtotals[i] + BABAtotals[i]; // assert(Ddenom3 == Ddenoms[i][2]);
        double D1 = Dnum1/Ddenom1; double D2 = Dnum2/Ddenom2; double D3 = Dnum3/Ddenom3;
        // Get the Z-scores
        double D1_Z = abs(D1)/D1stdErr; double D2_Z = abs(D2)/D2stdErr;
        double D3_Z = abs(D3)/D3stdErr;
        
        
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotals[i] >= BABAtotals[i] && BBAAtotals[i] >= ABBAtotals[i]) {
            if (D1 >= 0)
                *outFileBBAA << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2];
            else
                *outFileBBAA << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2];
            *outFileBBAA << "\t" << abs(D1) << "\t" << D1_Z << "\t";
            *outFileBBAA << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        } else if (BABAtotals[i] >= BBAAtotals[i] && BABAtotals[i] >= ABBAtotals[i]) {
            if (D2 >= 0)
                *outFileBBAA << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1];
            else
                *outFileBBAA << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1];
            *outFileBBAA << "\t" << abs(D2) << "\t" << D2_Z << "\t";
            *outFileBBAA << BABAtotals[i] << "\t" << BBAAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        } else if (ABBAtotals[i] >= BBAAtotals[i] && ABBAtotals[i] >= BABAtotals[i]) {
            if (D3 >= 0)
                *outFileBBAA << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0];
            else
                *outFileBBAA << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0];
            *outFileBBAA << "\t" << abs(D3) << "\t" << D3_Z << "\t";
            *outFileBBAA << ABBAtotals[i] << "\t" << BABAtotals[i] << "\t" << BBAAtotals[i] << std::endl;
        }
        
        // Find Dmin:
        if (abs(D1) <= abs(D2) && abs(D1) <= abs(D3)) { // (P3 == S3)
            if (D1 >= 0)
                *outFileDmin << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D1_Z << "\t" << std::endl;
            else
                *outFileDmin << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2] << "\t" << abs(D1) << "\t" << D1_Z << "\t"<< std::endl;
            // if (BBAAtotals[i] < BABAtotals[i] || BBAAtotals[i] < ABBAtotals[i])
            //     std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        } else if (abs(D2) <= abs(D1) && abs(D2) <= abs(D3)) { // (P3 == S2)
            if (D2 >= 0)
                *outFileDmin << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1] << "\t" << D2 << "\t" << D2_Z << "\t"<< std::endl;
            else
                *outFileDmin << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1] << "\t" << abs(D2) << "\t" << D2_Z << "\t"<< std::endl;
            // if (BABAtotals[i] < BBAAtotals[i] || BABAtotals[i] < ABBAtotals[i])
            //     std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        } else if (abs(D3) <= abs(D1) && abs(D3) <= abs(D2)) { // (P3 == S1)
            if (D3 >= 0)
                *outFileDmin << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0] << "\t" << D3 << "\t" << D3_Z << "\t"<< std::endl;
            else
                *outFileDmin << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0] << "\t" << abs(D3) << "\t" << D3_Z << "\t" << std::endl;;
            // if (ABBAtotals[i] < BBAAtotals[i] || ABBAtotals[i] < BABAtotals[i])
            //     std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        }
        
        // Output a simple file that can be used for combining multiple local runs:
        *outFileCombine << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        print_vector(regionDs[i][0], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t"; print_vector(regionDs[i][1], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t";
        print_vector(regionDs[i][2], *outFileCombineStdErr, ',',false); *outFileCombineStdErr << std::endl;
        
        //std::cerr << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D2 << "\t" << D3 << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
    }
     */
    return 0;
    
}



void parseDminCombineOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'h':
                std::cout << DMINCOMBINE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    
    int nFilenames = argc - optind;
    if (nFilenames < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DMINCOMBINE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    while (optind < argc) {
        opt::dminFiles.push_back(argv[optind++]);
    }
}
