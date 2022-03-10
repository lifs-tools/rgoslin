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

void assertEqual(string s1, string s2, string message = ""){
    if(s1 != s2){
        cout << "Assertion failed: '" << s1 << "' != '" << s2 << "'" << endl;
        if (message.length() > 0) cout << message << endl << endl;
        //exit(-1);
    }
}

int main(int argc, char** argv){
    string test_file = "data/goslin/testfiles/fatty-acids-test.csv";
    FattyAcidParser lipid_parser;
    ShorthandParser shorthand_parser;
        
    // test several more lipid names
    vector<string> lipid_data;
    vector<string> lipid_names;
    ifstream infile(test_file);
    assert(infile.good());
    string line;
    while (getline(infile, line)){
        line = strip(line, ' ');
        lipid_data.push_back(line);
    }
    infile.close();
    

    ////////////////////////////////////////////////////////////////////////////
    // Test for correctness
    ////////////////////////////////////////////////////////////////////////////
    
    
    for (auto lipid_row : lipid_data){
        
        vector<string> *data = split_string(lipid_row, ',', '"', true);
        string lmid = strip(data->at(0), '\"');
        string lipid_name = strip(data->at(1), '\"');
        string formula = strip(data->at(2), '\"');
        string expected_lipid_name = strip(data->at(3), '\"');
        
        LipidAdduct *lipid = 0;
        lipid = lipid_parser.parse(lipid_name);
        
        
        assertEqual(expected_lipid_name, lipid->get_lipid_string(), lmid + " '" + lipid_name + "': " + expected_lipid_name + " != " + lipid->get_lipid_string() + " (computed)");
            
        string lipid_formula = lipid->get_sum_formula();
        ElementTable *e = SumFormulaParser::get_instance().parse(formula);
        formula = goslin::compute_sum_formula(e);
        delete e;
        
        assertEqual(formula, lipid_formula, "formula " + lmid + " '" + lipid_name + "': " + formula + " != " + lipid_formula + " (computed)");
            
        if (to_lower(lipid_name).find("cyano") != string::npos) {
            delete lipid;
            delete data;
            continue;
        }
        
        LipidAdduct *lipid2 = shorthand_parser.parse(lipid->get_lipid_string());
        lipid_formula = lipid2->get_sum_formula();
        
        
        assertEqual(formula, lipid_formula, "lipid " + lmid + " '" + lipid_name + "': " + formula + " != " + lipid_formula + " (computed)");
        delete lipid2;
        
        lipid2 = shorthand_parser.parse(lipid->get_lipid_string(MOLECULAR_SPECIES));
        lipid_formula = lipid2->get_sum_formula();
        
        assertEqual(formula, lipid_formula, "molecular " + lmid + " '" + lipid_name + "': " + formula + " != " + lipid_formula + " (computed)");
        delete lipid2;
        
        lipid2 = shorthand_parser.parse(lipid->get_lipid_string(SPECIES));
        lipid_formula = lipid2->get_sum_formula();
        
        assertEqual(formula, lipid_formula, "species " + lmid + " '" + lipid_name + "': " + formula + " != " + lipid_formula + " (computed)");
        
        delete lipid2;
        delete lipid;
        delete data;
    }
    
    cout << "All tests passed without any problem" << endl;
    return 0;
}

