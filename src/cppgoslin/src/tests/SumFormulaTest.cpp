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
#include <cstdint>


using namespace std;
using namespace goslin;

int main(int argc, char** argv){
    
    LipidAdduct* lipid;
    SwissLipidsParser swiss_lipids_parser;
    LipidMapsParser lipid_maps_parser;
    
    /*
    lipid = swiss_lipids_parser.parse("SE(27:1/10:0)");
    cout << lipid->get_sum_formula() << " vs. C37H64O2" << endl;
    exit(0);
    */
    
    // test several more lipid names
    vector<string> lipid_names;
    vector<string> sum_formulas;
    ifstream infile("data/goslin/testfiles/formulas-lipid-maps.csv");
    assert(infile.good());
    string line;
    while (getline(infile, line)){
        vector<string>* tokens = split_string(line, ',', '"');
        line = strip(line, ' ');
        if (tokens->size() > 1){
            lipid_names.push_back(strip(tokens->at(0), '"'));
            sum_formulas.push_back(strip(tokens->at(1), '"'));
        }
        delete tokens;
    }
    infile.close();
        
    
    for (uint32_t i = 0; i < lipid_names.size(); ++i){
        string lipid_name = lipid_names.at(i);
        
        string correct_formula = sum_formulas.at(i);
        try {
            lipid = lipid_maps_parser.parse(lipid_name);
            assert(lipid != NULL);
            
            if (lipid->get_sum_formula() != correct_formula){
                cout << "Error for lipid '" << lipid_name << "': " << lipid->get_sum_formula() << " != " << correct_formula << " (reference)" << endl;
                assert(false);
            }
            delete lipid;
        }
        catch (LipidException &e){
            if (lipid_name.length() > 0){
                cout << "Exception: " << lipid_names.at(i) << endl;
                cout << e.what() << endl;
                assert(false);
            }
        }
    }
    cout << "All tests passed without any problem" << endl;
    
    
    
    // test several more lipid names
    lipid_names.clear();
    sum_formulas.clear();
    
    ifstream infile2("data/goslin/testfiles/formulas-swiss-lipids.csv");
    assert(infile2.good());
    while (getline(infile2, line)){
        vector<string>* tokens = split_string(line, ',', '"');
        line = strip(line, ' ');
        if (tokens->size() > 1){
            lipid_names.push_back(strip(tokens->at(0), '"'));
            sum_formulas.push_back(strip(tokens->at(1), '"'));
        }
        delete tokens;
    }
    infile2.close();
        
    
    for (uint32_t i = 0; i < lipid_names.size(); ++i){
        string lipid_name = lipid_names.at(i);
        string correct_formula = sum_formulas.at(i);
        //cout << lipid_name << " " << correct_formula << endl;
        try {
            lipid = swiss_lipids_parser.parse(lipid_name);
            assert(lipid != NULL);
            
            if (lipid->get_sum_formula() != correct_formula){
                cout << "Error for lipid '" << lipid_name << "': " << lipid->get_sum_formula() << " != " << correct_formula << endl;
                assert(false);
            }
            delete lipid;
        }
        catch (LipidException &e){
            if (lipid_name.length() > 0){
                cout << "Exception: " << lipid_names.at(i) << endl;
                cout << e.what() << endl;
                assert(false);
            }
        }
    }
    
    cout << "All tests passed without any problem" << endl;

    return EXIT_SUCCESS;
}
