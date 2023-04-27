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
        cout << "Assertion: '" << a << "' == '" << b << "'" << endl;
        assert(a == b);
    }
}

void assert_true(int a, int b, string t = "") {
    if (a != b){
        cout << "Assertion: '" << a << "' == '" << b << "'" << endl;
        assert(a == b);
    }
}

void assert_true(std::size_t a, std::size_t b, string t = "") {
    if (a != b){
        cout << "Assertion: '" << a << "' == '" << b << "'" << endl;
        assert(a == b);
    }
}

int main(int argc, char** argv){
    
    LipidAdduct* lipid;
    string test_file = "data/goslin/testfiles/goslin-test.csv";
    GoslinParser parser;
        
    LipidAdduct *l = parser.parse("Cer 18:1(8Z);2/24:0");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "Cer 18:1(8);(OH)2/24:0");
    assert_true(l->get_lipid_string(SN_POSITION), "Cer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "Cer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(SPECIES), "Cer 42:1;O2");
    assert_true(l->get_sum_formula(), "C42H83NO3");
    
    l = parser.parse("HexCer 18:1(5Z);2/24:0");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "HexCer 18:1(5);OH/24:0");
    assert_true(l->get_lipid_string(SN_POSITION), "HexCer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "HexCer 18:1;O2/24:0");
    assert_true(l->get_lipid_string(SPECIES), "HexCer 42:1;O2");
    assert_true(l->get_sum_formula(), "C48H93NO8");
    
    l = parser.parse("LSM 17:1(4E);2");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "LSM 17:1(4);OH");
    assert_true(l->get_lipid_string(SN_POSITION), "LSM 17:1;O2");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "LSM 17:1;O2");
    assert_true(l->get_lipid_string(SPECIES), "LSM 17:1;O2");
    assert_true(l->get_sum_formula(), "C22H47N2O5P");
    
    l = parser.parse("LCB 18:1(4E);2");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "SPB 18:1(4);(OH)2");
    assert_true(l->get_lipid_string(SN_POSITION), "SPB 18:1;O2");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "SPB 18:1;O2");
    assert_true(l->get_lipid_string(SPECIES), "SPB 18:1;O2");
    assert_true(l->get_sum_formula(), "C18H37NO2");

    l = parser.parse("EPC 14:1(4E);2/20:1(11Z)");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "EPC 14:1(4);OH/20:1(11)");
    assert_true(l->get_lipid_string(SN_POSITION), "EPC 14:1;O2/20:1");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "EPC 14:1;O2/20:1");
    assert_true(l->get_lipid_string(SPECIES), "EPC 34:2;O2");
    assert_true(l->get_sum_formula(), "C36H71N2O6P");
    
    l = parser.parse("MIPC 18:0;3/24:0");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "MIPC 18:0;(OH)2/24:0");
    assert_true(l->get_lipid_string(SN_POSITION), "MIPC 18:0;O3/24:0");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "MIPC 18:0;O3/24:0");
    assert_true(l->get_lipid_string(SPECIES), "MIPC 42:0;O3");
    assert_true(l->get_sum_formula(), "C54H106NO17P");
    
    l = parser.parse("EPC 16:2(4E,6E);2/22:1(13Z);1");
    assert_true(l->get_lipid_string(STRUCTURE_DEFINED), "EPC 16:2(4,6);OH/22:1(13);OH");
    assert_true(l->get_lipid_string(SN_POSITION), "EPC 16:2;O2/22:1;O");
    assert_true(l->get_lipid_string(MOLECULAR_SPECIES), "EPC 16:2;O2/22:1;O");
    assert_true(l->get_lipid_string(SPECIES), "EPC 38:3;O3");
    assert_true(l->get_sum_formula(), "C40H77N2O7P");

    l = parser.parse("BMP 18:1-18:1");
    assert_true(l->get_lipid_string(), "BMP 18:1_18:1");
    assert_true(l->get_sum_formula(), "C42H79O10P");
    assert_true(l->lipid->fa_list.at(0)->name, "FA1");
    assert_true(l->lipid->fa_list.at(0)->position, -1);
    assert_true(l->lipid->fa_list.at(0)->double_bonds->num_double_bonds, 1);
    assert_true(l->lipid->fa_list.at(1)->name, "FA2");
    assert_true(l->lipid->fa_list.at(1)->position, -1);
    assert_true(l->lipid->fa_list.at(1)->double_bonds->num_double_bonds, 1);
    assert_true(l->lipid->fa_list.at(2)->name,  "FA3");
    assert_true(l->lipid->fa_list.at(2)->position, -1);
    assert_true(l->lipid->fa_list.at(2)->num_carbon, 0);
    assert_true(l->lipid->fa_list.at(2)->double_bonds->num_double_bonds, 0);
    assert_true(l->lipid->fa_list.at(3)->name,  "FA4");
    assert_true(l->lipid->fa_list.at(3)->position, -1);
    assert_true(l->lipid->fa_list.at(3)->num_carbon, 0);
    assert_true(l->lipid->fa_list.at(3)->double_bonds->num_double_bonds, 0);

    l = parser.parse("TAG 18:1/0:0/16:0");
    assert_true(l->get_lipid_string(), "DG 18:1/0:0/16:0");
    assert_true(l->get_sum_formula(), "C37H70O5");
    assert_true(l->lipid->fa_list.size(), (std::size_t)3);
    assert_true(l->lipid->fa_list.at(0)->name, "FA1");
    assert_true(l->lipid->fa_list.at(0)->position, 1);
    assert_true(l->lipid->fa_list.at(0)->num_carbon, 18);
    assert_true(l->lipid->fa_list.at(0)->double_bonds->num_double_bonds, 1);
    assert_true(l->lipid->fa_list.at(1)->name, "FA2");
    assert_true(l->lipid->fa_list.at(1)->position, 2);
    assert_true(l->lipid->fa_list.at(1)->num_carbon, 0);
    assert_true(l->lipid->fa_list.at(1)->double_bonds->num_double_bonds, 0);
    assert_true(l->lipid->fa_list.at(2)->name, "FA3");
    assert_true(l->lipid->fa_list.at(2)->position, 3);
    assert_true(l->lipid->fa_list.at(2)->num_carbon, 16);
    assert_true(l->lipid->fa_list.at(2)->double_bonds->num_double_bonds, 0);
    
    l = parser.parse("15S-HETE-d8");
    assert_true(l->get_lipid_string(), "FA 20:4;OH[M[2]H8]");
    assert_true(l->get_sum_formula(), "C20H24O3H'8");
    
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
        try {
            lipid = parser.parse(lipid_name);
            assert(lipid != NULL);
            delete lipid;
        }
        catch (LipidException &e){
            cout << "Exception: " << lipid_name << endl;
            cout << e.what() << endl;
        }
    }
    
    
    cout << "All tests passed without any problem" << endl;

    return EXIT_SUCCESS;
}
