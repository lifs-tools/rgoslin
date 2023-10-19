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
#include <fstream>
#include <vector>
#include <map>
#include <cassert>


using namespace std;
using namespace goslin;


       
        
void assertEqual(string s1, string s2, string message = ""){
    if(s1 != s2){
        cout << "Assertion failed: '" << s1 << "' != '" << s2 << "' (reference)" << endl;
        if (message.length() > 0) cout << message << endl;
        exit(-1);
    }
}
        

int main(int argc, char** argv){
    ShorthandParser parser;
        
    LipidAdduct *l = parser.parse("PE 18:1(8Z);1OH,3OH/24:0");
    assertEqual(l->get_lipid_string(), "PE 18:1(8Z);1OH,3OH/24:0");
    assertEqual(l->get_lipid_string(STRUCTURE_DEFINED), "PE 18:1(8);(OH)2/24:0");
    assertEqual(l->get_lipid_string(SN_POSITION), "PE 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "PE 18:1;O2_24:0");
    assertEqual(l->get_lipid_string(SPECIES), "PE 42:1;O2");
    delete l;
        
    l = parser.parse("Cer 18:1(8Z);1OH,3OH/24:0");
    assertEqual(l->get_lipid_string(), "Cer 18:1(8Z);1OH,3OH/24:0");
    assertEqual(l->get_lipid_string(STRUCTURE_DEFINED), "Cer 18:1(8);(OH)2/24:0");
    assertEqual(l->get_lipid_string(SN_POSITION), "Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(SPECIES), "Cer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C42H83NO3");
    delete l;
    
    l = parser.parse("Cer 18:1(8);(OH)2/24:0");
    assertEqual(l->get_lipid_string(), "Cer 18:1(8);(OH)2/24:0");
    assertEqual(l->get_lipid_string(SN_POSITION), "Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(SPECIES), "Cer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C42H83NO3");
    delete l;
    
    l = parser.parse("Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(), "Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "Cer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(SPECIES), "Cer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C42H83NO3");
    delete l;
    
    l = parser.parse("Cer 42:1;O2");
    assertEqual(l->get_lipid_string(), "Cer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C42H83NO3");
    delete l;
    
    l = parser.parse("Gal-Cer(1) 18:1(5Z);3OH/24:0");
    assertEqual(l->get_lipid_string(), "Gal-Cer(1) 18:1(5Z);3OH/24:0");
    assertEqual(l->get_lipid_string(STRUCTURE_DEFINED), "Gal-Cer 18:1(5);OH/24:0");
    assertEqual(l->get_lipid_string(SN_POSITION), "HexCer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "HexCer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(SPECIES), "HexCer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C48H93NO8");
    delete l;
    
    l = parser.parse("Gal-Cer 18:1(5);OH/24:0");
    assertEqual(l->get_lipid_string(), "Gal-Cer 18:1(5);OH/24:0");
    assertEqual(l->get_lipid_string(SN_POSITION), "HexCer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "HexCer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(SPECIES), "HexCer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C48H93NO8");
    delete l;
    
    
    l = parser.parse("HexCer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(), "HexCer 18:1;O2/24:0");
    assertEqual(l->get_lipid_string(SPECIES), "HexCer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C48H93NO8");
    delete l;
    
    
    l = parser.parse("HexCer 42:1;O2");
    assertEqual(l->get_lipid_string(), "HexCer 42:1;O2");
    assertEqual(l->get_sum_formula(), "C48H93NO8");
    delete l;
    
    
    
    l = parser.parse("SPB 18:1(4Z);1OH,3OH");
    assertEqual(l->get_lipid_string(), "SPB 18:1(4Z);1OH,3OH");
    assertEqual(l->get_lipid_string(STRUCTURE_DEFINED), "SPB 18:1(4);(OH)2");
    assertEqual(l->get_lipid_string(SN_POSITION), "SPB 18:1;O2");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "SPB 18:1;O2");
    assertEqual(l->get_lipid_string(SPECIES), "SPB 18:1;O2");
    assertEqual(l->get_sum_formula(), "C18H37NO2");
    delete l;
    
    l = parser.parse("SPB 18:1(4);(OH)2");
    assertEqual(l->get_lipid_string(), "SPB 18:1(4);(OH)2");
    assertEqual(l->get_lipid_string(SN_POSITION), "SPB 18:1;O2");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "SPB 18:1;O2");
    assertEqual(l->get_lipid_string(SPECIES), "SPB 18:1;O2");
    assertEqual(l->get_sum_formula(), "C18H37NO2");
    delete l;
    
    l = parser.parse("SPB 18:1;O2");
    assertEqual(l->get_lipid_string(), "SPB 18:1;O2");
    assertEqual(l->get_lipid_string(SPECIES), "SPB 18:1;O2");
    assertEqual(l->get_sum_formula(), "C18H37NO2");
    delete l;
    

    
    l = parser.parse("LSM(1) 17:1(4E);3OH");
    assertEqual(l->get_lipid_string(), "LSM(1) 17:1(4E);3OH");
    assertEqual(l->get_lipid_string(STRUCTURE_DEFINED), "LSM 17:1(4);OH");
    assertEqual(l->get_lipid_string(SN_POSITION), "LSM 17:1;O2");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "LSM 17:1;O2");
    assertEqual(l->get_lipid_string(SPECIES), "LSM 17:1;O2");
    assertEqual(l->get_sum_formula(), "C22H47N2O5P");
    delete l;
    
    l = parser.parse("LSM 17:1(4);OH");
    assertEqual(l->get_lipid_string(), "LSM 17:1(4);OH");
    assertEqual(l->get_lipid_string(SN_POSITION), "LSM 17:1;O2");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "LSM 17:1;O2");
    assertEqual(l->get_lipid_string(SPECIES), "LSM 17:1;O2");
    assertEqual(l->get_sum_formula(), "C22H47N2O5P");
    delete l;
    
    l = parser.parse("LSM 17:1;O2");
    assertEqual(l->get_lipid_string(), "LSM 17:1;O2");
    assertEqual(l->get_lipid_string(SPECIES), "LSM 17:1;O2");
    assertEqual(l->get_sum_formula(), "C22H47N2O5P");
    delete l;
    

    
    l = parser.parse("EPC(1) 14:1(4E);3OH/20:1(11Z)");
    assertEqual(l->get_lipid_string(), "EPC(1) 14:1(4E);3OH/20:1(11Z)");
    assertEqual(l->get_lipid_string(STRUCTURE_DEFINED), "EPC 14:1(4);OH/20:1(11)");
    assertEqual(l->get_lipid_string(SN_POSITION), "EPC 14:1;O2/20:1");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "EPC 14:1;O2/20:1");
    assertEqual(l->get_lipid_string(SPECIES), "EPC 34:2;O2");
    assertEqual(l->get_sum_formula(), "C36H71N2O6P");
    delete l;
    
    l = parser.parse("EPC 14:1(4);OH/20:1(11)");
    assertEqual(l->get_lipid_string(), "EPC 14:1(4);OH/20:1(11)");
    assertEqual(l->get_lipid_string(SN_POSITION), "EPC 14:1;O2/20:1");
    assertEqual(l->get_lipid_string(MOLECULAR_SPECIES), "EPC 14:1;O2/20:1");
    assertEqual(l->get_lipid_string(SPECIES), "EPC 34:2;O2");
    assertEqual(l->get_sum_formula(), "C36H71N2O6P");
    delete l;
    
    l = parser.parse("EPC 14:1;O2/20:1");
    assertEqual(l->get_lipid_string(), "EPC 14:1;O2/20:1");
    assertEqual(l->get_lipid_string(SPECIES), "EPC 34:2;O2");
    assertEqual(l->get_sum_formula(), "C36H71N2O6P");
    delete l;
    
    l = parser.parse("EPC 34:2;O2");
    assertEqual(l->get_lipid_string(), "EPC 34:2;O2");
    assertEqual(l->get_sum_formula(), "C36H71N2O6P");
    delete l;
    
    l = parser.parse("SPBP 18:1;O2");
    assertEqual(l->get_lipid_string(), "SPBP 18:1;O2");
    assertEqual(l->get_sum_formula(), "C18H38NO5P");
    delete l;
        
    l = parser.parse("LSM 18:1;O2");
    assertEqual(l->get_lipid_string(), "LSM 18:1;O2");
    delete l;
    
    l = parser.parse("LSM(1) 18:1(5Z);3OH");
    assertEqual(l->get_lipid_string(), "LSM(1) 18:1(5Z);3OH");
    delete l;
    
    l = parser.parse("LHexCer 18:1;O2");
    assertEqual(l->get_lipid_string(), "LHexCer 18:1;O2");
    delete l;
    
    l = parser.parse("LHexCer 18:1;O2/0:0");
    assertEqual(l->get_lipid_string(), "LHexCer 18:1;O2");
    delete l;
    
    l = parser.parse("LHexCer(1) 18:1(5E);3OH/0:0");
    assertEqual(l->get_lipid_string(), "LHexCer(1) 18:1(5E);3OH");
    delete l;
    
    l = parser.parse("TG 42:2_18:0");
    assertEqual(l->get_lipid_string(), "TG 42:2_18:0");
    delete l;
    
    
    // test several more lipid names
    vector<string> data;
    vector<LipidLevel> levels {FULL_STRUCTURE, STRUCTURE_DEFINED, SN_POSITION, MOLECULAR_SPECIES, SPECIES};
    
    ifstream infile("data/goslin/testfiles/shorthand-test.csv");
    assert(infile.good());
    string line;
    while (getline(infile, line)) data.push_back(line);
    infile.close();
    
    try {
        l = parser.parse("SM 21:1(3Z);2O/12:0");
        assert(false);
    }
    catch (LipidException &){
        // great, we wanted the previous name to fail
    }
    
    
    for (auto &row : data){
        vector<string>* results = split_string(row, ',', '"', true);
        
        //cout << results->at(0) << endl;
        for (int i = 0; i < (int)results->size(); ++i) results->at(i) = strip(results->at(i), '"');
        string lipid_name = results->at(0);
        
        LipidAdduct *lipid = parser.parse(lipid_name);
        string formula = (results->size() > levels.size() && results->at(levels.size()).length() > 0) ? results->at(levels.size()) : lipid->get_sum_formula();
        
        
        if (results->size() > levels.size()){
            assertEqual(formula, lipid->get_sum_formula(), "test on lipid '" + lipid_name + "'");
        }
        
        for (int l = 0; l < (int)levels.size(); ++l){
            LipidLevel lipid_level = levels.at(l);
            string n = lipid->get_lipid_string(lipid_level);
            assertEqual(results->at(l), n, "test 1 on lipid '" + lipid_name + "' and level '" + std::to_string(lipid_level) + "'");
            assertEqual(formula, lipid->get_sum_formula(), "test 2 on lipid '" + lipid_name + "' and level '" + std::to_string(lipid_level) + "'");

            LipidAdduct *lipid2 = parser.parse(n);
            for (int ll = l; ll < (int)levels.size(); ++ll){
                assertEqual(results->at(ll), lipid2->get_lipid_string(levels.at(ll)), "test 3 on lipid '" + n + "' / " + results->at(ll) + " vs " + lipid2->get_lipid_string(levels.at(ll)) + " and level '" + std::to_string(levels.at(ll)) + "'");
                assertEqual(formula, lipid2->get_sum_formula(), "test 4 on lipid '" + lipid_name + "' and level '" + std::to_string(levels.at(l)) + "' mapped to level '" + std::to_string(levels.at(ll)) + "'");
            }
            delete lipid2;
        }
        delete lipid;
        delete results;
    }
    cout << "All tests passed without any problem" << endl;
    return 0;
}
