/*
MIT License

Copyright (c) the authors (listed in global LICENSE file)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include "cppgoslin/cppgoslin.h"
#include <set>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
#include <map>


using namespace std;
using namespace goslin;


void assert_true(string a, string b, string t = ""){
    if (a != b){
        cout << "Assertion: '" << a << "' == '" << b << " (" << t << ")" << endl;
        assert(a == b);
    }
}


int main(int argc, char** argv){
    
    LipidAdduct* lipid;
    LipidMapsParser parser;
    
    
    
    lipid = parser.parse("Cer(d18:1(8Z)/24:0)");
    assert_true(lipid->get_lipid_string(STRUCTURE_DEFINED), "Cer 18:1(8);(OH)2/24:0");
    assert_true(lipid->get_lipid_string(SN_POSITION), "Cer 18:1;O2/24:0");
    assert_true(lipid->get_lipid_string(MOLECULAR_SPECIES), "Cer 18:1;O2/24:0");
    assert_true(lipid->get_lipid_string(SPECIES), "Cer 42:1;O2");
    assert_true(lipid->get_sum_formula(), "C42H83NO3");
    delete lipid;
    
    lipid = parser.parse("GalCer(d18:1(5Z)/24:0)");
    assert_true(lipid->get_lipid_string(STRUCTURE_DEFINED), "GalCer 18:1(5);OH/24:0");
    assert_true(lipid->get_lipid_string(SN_POSITION), "GalCer 18:1;O2/24:0");
    assert_true(lipid->get_lipid_string(MOLECULAR_SPECIES), "GalCer 18:1;O2/24:0");
    assert_true(lipid->get_lipid_string(SPECIES), "GalCer 42:1;O2");
    assert_true(lipid->get_sum_formula(), "C48H93NO8");
    delete lipid;
    
    lipid = parser.parse("LysoSM(d17:1(4E))");
    assert_true(lipid->get_lipid_string(STRUCTURE_DEFINED), "LSM 17:1(4);OH");
    assert_true(lipid->get_lipid_string(SN_POSITION), "LSM 17:1;O2");
    assert_true(lipid->get_lipid_string(MOLECULAR_SPECIES), "LSM 17:1;O2");
    assert_true(lipid->get_lipid_string(SPECIES), "LSM 17:1;O2");
    assert_true(lipid->get_sum_formula(), "C22H47N2O5P");
    delete lipid;

    lipid = parser.parse("PE-Cer(d14:1(4E)/20:1(11Z))");
    assert_true(lipid->get_lipid_string(STRUCTURE_DEFINED), "EPC 14:1(4);OH/20:1(11)");
    assert_true(lipid->get_lipid_string(SN_POSITION), "EPC 14:1;O2/20:1");
    assert_true(lipid->get_lipid_string(MOLECULAR_SPECIES), "EPC 14:1;O2/20:1");
    assert_true(lipid->get_lipid_string(SPECIES), "EPC 34:2;O2");
    assert_true(lipid->get_sum_formula(), "C36H71N2O6P");
    delete lipid;
    
    lipid = parser.parse("MIPC(t18:0/24:0)");
    assert_true(lipid->get_lipid_string(STRUCTURE_DEFINED), "MIPC 18:0;(OH)2/24:0");
    assert_true(lipid->get_lipid_string(SN_POSITION), "MIPC 18:0;O3/24:0");
    assert_true(lipid->get_lipid_string(MOLECULAR_SPECIES), "MIPC 18:0;O3/24:0");
    assert_true(lipid->get_lipid_string(SPECIES), "MIPC 42:0;O3");
    assert_true(lipid->get_sum_formula(), "C54H106NO17P");
    delete lipid;
    
    lipid = parser.parse("PE-Cer(d16:2(4E,6E)/22:1(13Z)(2OH))");
    assert_true(lipid->get_lipid_string(STRUCTURE_DEFINED), "EPC 16:2(4,6);OH/22:1(13);OH");
    assert_true(lipid->get_lipid_string(SN_POSITION), "EPC 16:2;O2/22:1;O");
    assert_true(lipid->get_lipid_string(MOLECULAR_SPECIES), "EPC 16:2;O2/22:1;O");
    assert_true(lipid->get_lipid_string(SPECIES), "EPC 38:3;O3");
    assert_true(lipid->get_sum_formula(), "C40H77N2O7P");
    delete lipid;

    lipid = parser.parse("GalNAcβ1-4(Galβ1-4GlcNAcβ1-3)Galβ1-4Glcβ-Cer(d18:1/24:1(15Z))");
    assert_true(lipid->get_lipid_string(SPECIES), "GalNAcGalGlcNAcGalGlcCer 42:2;O2");
    delete lipid;
    
    // test several more lipid names
    vector<string> lipid_names_income;
    vector<string> lipid_names_outcome;
    ifstream infile("data/goslin/testfiles/lipid-maps-test.csv");
    assert(infile.good());
    string line;
    while (getline(infile, line)){
        vector<string>* tokens = split_string(line, ',', '"', true);
        line = strip(line, ' ');
        if (tokens->size() > 1){
            lipid_names_income.push_back(strip(tokens->at(0), '"'));
            lipid_names_outcome.push_back(strip(tokens->at(1), '"'));
        }
        delete tokens;
    }
    infile.close();
        
    #pragma omp parallel for
    for (int i = 0; i < (int)lipid_names_income.size(); ++i){
        string lipid_name = lipid_names_income.at(i);
        string correct_lipid_name = lipid_names_outcome.at(i);
        //cout << lipid_name << endl;
        try {
            LipidAdduct* l = parser.parse_parallel(lipid_name);
            assert(l != NULL);
            if (correct_lipid_name != "Unsupported lipid" && correct_lipid_name.length() > 0){
                assert_true(correct_lipid_name, l->get_lipid_string(), lipid_name + ": " + l->get_lipid_string() + " != " + correct_lipid_name + " (reference)");
            }
            delete l;
        }
        catch (exception &e){
            if (correct_lipid_name.length() > 0){
                cout << "Exception: " << i << ": " << lipid_names_income.at(i) << endl;
                cout << e.what() << endl;
                assert(false);
            }
        }
    }
    
    
    cout << "All tests passed without any problem" << endl;

    
    
    return EXIT_SUCCESS;
}
