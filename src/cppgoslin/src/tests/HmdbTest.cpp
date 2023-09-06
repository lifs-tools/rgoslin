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


using namespace std;
using namespace goslin;


void assert_true(string a, string b, string t = ""){
    if (a != b){
        cout << "Assertion: '" << a << "' == '" << b << endl;
        assert(a == b);
    }
}

int main(int argc, char** argv){
    
    LipidAdduct* lipid;
    string test_file = "data/goslin/testfiles/hmdb-test.csv";
    HmdbParser parser;
    
    
    LipidAdduct* l = parser.parse("Cer(d18:1(8Z)/24:0)");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "Cer 18:1(8);(OH)2/24:0");
    assert_true(l->get_lipid_string(SN_POSITION), "Cer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "Cer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(SPECIES), "Cer 42:1;O2");
    assert_true(l->get_sum_formula(), "C42H83NO3");
    delete l;
    
    l = parser.parse("GalCer(d18:1(5Z)/24:0)");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "GalCer 18:1(5);OH/24:0");
    assert_true(l->get_lipid_string(SN_POSITION), "GalCer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "GalCer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(SPECIES), "GalCer 42:1;O2");
    assert_true(l->get_sum_formula(), "C48H93NO8");
    delete l;
    
    l = parser.parse("LysoSM(d17:1(4E))");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "LSM 17:1(4);OH");
    assert_true(l->get_lipid_string(SN_POSITION), "LSM 17:1;O2");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "LSM 17:1;O2");
    assert_true(l->get_lipid_string(SPECIES), "LSM 17:1;O2");
    assert_true(l->get_sum_formula(), "C22H47N2O5P");
    delete l;

    l = parser.parse("PE-Cer(d14:1(4E)/20:1(11Z))");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "EPC 14:1(4);OH/20:1(11)");
    assert_true(l->get_lipid_string(SN_POSITION), "EPC 14:1;O2/20:1");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "EPC 14:1;O2/20:1");
    assert_true(l->get_lipid_string(SPECIES), "EPC 34:2;O2");
    assert_true(l->get_sum_formula(), "C36H71N2O6P");
    delete l;
    
    l = parser.parse("MIPC(t18:0/24:0)");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "MIPC 18:0;(OH)2/24:0");
    assert_true(l->get_lipid_string(SN_POSITION), "MIPC 18:0;O3/24:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "MIPC 18:0;O3/24:0");
    assert_true(l->get_lipid_string(SPECIES), "MIPC 42:0;O3");
    assert_true(l->get_sum_formula(), "C54H106NO17P");
    delete l;
    
    l = parser.parse("PE-Cer(d16:2(4E,6E)/22:1(13Z)(2OH))");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "EPC 16:2(4,6);OH/22:1(13);OH");
    assert_true(l->get_lipid_string(SN_POSITION), "EPC 16:2;O2/22:1;O");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "EPC 16:2;O2/22:1;O");
    assert_true(l->get_lipid_string(SPECIES), "EPC 38:3;O3");
    assert_true(l->get_sum_formula(), "C40H77N2O7P");
    delete l;
    
    l = parser.parse("DG(a-21:0/20:5(5Z,8Z,10E,14Z,17Z)+=O(12S)/0:0)");
    assert_true(l->get_lipid_string(FULL_STRUCTURE), "DG 20:0;Me/20:5(5Z,8Z,10E,14Z,17Z);12oxo/0:0");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "DG 20:0;Me/20:5(5,8,10,14,17);oxo/0:0");
    assert_true(l->get_lipid_string(SN_POSITION), "DG 21:0/20:6;O/0:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "DG 21:0_20:6;O");
    assert_true(l->get_lipid_string(SPECIES), "DG 41:6;O");
    assert_true(l->get_sum_formula(), "C44H74O6");
    delete l;
    
    l = parser.parse("PIP(16:2(9Z,12Z)/PGJ2)");
    assert_true(l->get_lipid_string(SN_POSITION), "PIP 16:2/20:5;O2");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "PIP 16:2_20:5;O2");
    assert_true(l->get_lipid_string(SPECIES), "PIP 36:7;O2");
    assert_true(l->get_sum_formula(), "C45H74O18P2");
    delete l;
    
    l = parser.parse("PC(14:0/20:5(7Z,9Z,11E,13E,17Z)-3OH(5,6,15))");
    assert_true(l->get_lipid_string(FULL_STRUCTURE), "PC 14:0/20:5(7Z,9Z,11E,13E,17Z);5OH,6OH,15OH");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "PC 14:0/20:5(7,9,11,13,17);(OH)3");
    assert_true(l->get_lipid_string(SN_POSITION), "PC 14:0/20:5;O3");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "PC 14:0_20:5;O3");
    assert_true(l->get_lipid_string(SPECIES), "PC 34:5;O3");
    assert_true(l->get_sum_formula(), "C42H74NO11P");
    delete l;
    
    
        
    // test several more lipid names
    vector<string> lipid_names;
    ifstream infile(test_file);
    assert(infile.good());
    string line;
    while (getline(infile, line)){
        line = strip(line, ' ');
        lipid_names.push_back(line);
    }
    infile.close();
    
    for (auto lipid_name : lipid_names){
        lipid = parser.parse(lipid_name);
        assert(lipid != 0);
    }
    
    cout << "All tests passed without any problem" << endl;
    
    return EXIT_SUCCESS;
}
