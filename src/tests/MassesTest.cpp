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
#include <math.h>
#include <map>
#include <cassert>


using namespace std;
using namespace goslin;


void cpp_assert(string s1, string s2, string comment = ""){
    if (s1 != s2){
        cout << "assertion failed: " << s1 << " != " << s2 << " " << comment << endl;
        assert(false);
    }
}

void cpp_assert(double d1, double d2){
    if (fabs(d1 - d2) > 1e-4){
        cout << "assertion failed: " << d1 << " != " << d2 << endl;
        assert(false);
    }
}




int main(int argc, char** argv){
    
    LipidAdduct* lipid;
    GoslinParser parser;
    
    
    // test several more lipid names
    ifstream infile("data/goslin/testfiles/lipid-masses.csv");
    assert(infile.good());
    string line;
    int ii = 0;
    while (getline(infile, line)){
        if (ii++ == 0) continue;
        
        vector<string>* tokens = split_string(line, ',', '"');
        line = strip(line, ' ');
        if (tokens->size() > 1){
            string lipid_class = strip(tokens->at(0), '"');
            string lipid_name = strip(tokens->at(1), '"');
            string lipid_formula = strip(tokens->at(2), '"');
            string lipid_adduct = strip(tokens->at(3), '"');
            double lipid_mass = atof(strip(tokens->at(4), '"').c_str());
            int lipid_charge = atoi(strip(tokens->at(5), '"').c_str());
            
            string full_lipid_name = lipid_name + lipid_adduct;
            try{
                lipid = parser.parse(full_lipid_name);
                assert(lipid != NULL);
            }
            catch (LipidException& e){
                cout << lipid_name << ": " << e.what() << endl;
                assert (false);
            }
            
            cpp_assert (lipid->get_lipid_string(CLASS), lipid_class, " lipid name: " + lipid_name);
            cpp_assert (compute_sum_formula(lipid->lipid->get_elements()), lipid_formula);
            cpp_assert (lipid->get_mass(), lipid_mass);
            cpp_assert (lipid->adduct->get_charge(), lipid_charge);
            
            delete lipid;
            
        }
        delete tokens;
    }
    infile.close();
    
    
    cout << "All tests passed without any problem" << endl;

    
    
    return EXIT_SUCCESS;
}
