/*
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

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


#include <iostream>
#include "cppgoslin/cppgoslin.h"
#include <string>
#include <vector>
#include <map>

using namespace std;
using namespace goslin;

int main(){
    try {
        /* initiating the parser class */
        LipidParser lipid_parser;
        LipidAdduct* lipid;
            
        /* read in file with lipid names, one per row */
        vector<string> lipidnames;
        ifstream infile("../data/goslin/testfiles/lipidnames.txt");
        if (!infile.is_open()){
            cout << "Cannot find example file 'lipidnames.txt', exit." << endl;
            return -1;
        }
        string line;
        while (getline(infile, line)){
            line = strip(line, ' ');
            if (line.length() < 2) continue;
            lipidnames.push_back(line);
        }
        infile.close();
        
        /* setting up a map for counting the distribution of
         * the lipids within the lipid categories
         */
        map<string, int> lipid_counts;
        
        /* parsing all lipids */
        for (auto lipid_name : lipidnames){
            lipid = lipid_parser.parse(lipid_name);
            
            /* checking if lipid name was parsed */
            if (lipid){
                string category = lipid->get_lipid_string(CATEGORY);
                delete lipid;
                
                /* adding category into the map */
                if (lipid_counts.find(category) == lipid_counts.end()) lipid_counts.insert({category, 0});
                ++lipid_counts.at(category);
            }
        }
        
        /* reporting the distribution */
        cout << "Lipid distribution:" << endl;
        for (auto key_value_pair : lipid_counts){
            cout << key_value_pair.first << ": " << key_value_pair.second << " lipids" << endl;
        }
    }
    catch (LipidException &e){
        cout << "Exception:" << endl;
        cout << e.what() << endl;
    }
    
    
    return 0;
}
