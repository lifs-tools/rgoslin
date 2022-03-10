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
#include <fstream>

using namespace std;
using namespace goslin;

int main(){
  try {
    /* initiating the parser class */
    SwissLipidsParser lipid_parser;
    LipidAdduct* lipidAdduct = NULL;
    
    /* read in file with lipid names, one per row */
    vector<string> lipidnames;
    ifstream infile("../data/goslin/examplelists/swisslipids-names-only.tsv");
    if (!infile.is_open()){
        cout << "Cannot find example file 'swisslipids-names-only.tsv', exit." << endl;
        return -1;
    }
    string line;
    while (getline(infile, line)){
      line = strip(line, ' ');
      if (line.length() < 2) continue;
      lipidnames.push_back(line);
    }
    infile.close();
    cout << "Reading input complete" << endl;  
    /* setting up a map for counting the distribution of
     * the lipids within the lipid categories
     */
    map<string, int> lipid_counts;
    
    cout << "Loaded " << lipidnames.size() << " lipid names." << endl;
    
    cout << "Parsing lipids" << endl; 
    int parsedLipids = 0;    
    int totalLipids = 0;
    ofstream slout;
    string fileName = "swisslipids-results.tsv";
    slout.open(fileName, std::ofstream::out | std::ofstream::trunc);
    if(slout.is_open()) {
      cout << "Saving parsing results to '" << fileName << "'" << endl; 
      slout << "Original Name" << '\t' << "Structural Level" << '\t' << "Name at Level" << '\t' << "Category" << '\t' << "Species" << endl;
      /* parsing all lipids */
      for (auto lipid_name : lipidnames){
        ++totalLipids;
        slout << lipid_name << '\t';
        try {
          cout << "Parsing lipid " << lipid_name << endl;
          lipidAdduct = lipid_parser.parse(lipid_name);
        } catch (LipidException &e){
          cout << "Exception while parsing lipid " << lipid_name << ":" << endl;
          cout << e.what() << endl;
        }
        /* checking if lipid name was parsed */
        if (lipidAdduct != NULL){
          ++parsedLipids;
          LipidSpeciesInfo *info = lipidAdduct->lipid->info;
          string nativeLevelName = "";
          string category = "";
          string species = "";
          try {
            nativeLevelName = lipidAdduct->get_lipid_string(info->level);
          } catch (ConstraintViolationException &ce) {
            cout << "Constraint violation:" << endl;
            cout << ce.what() << endl;
          } catch (RuntimeException &re) {
            cout << "Runtime exception while generating native level name for lipid " << lipid_name << endl;
          }
          slout << info->level << '\t';
          slout << nativeLevelName << '\t';
          try {            
            category = lipidAdduct->get_lipid_string(CATEGORY);
          } catch (ConstraintViolationException &ce) {
            cout << "Constraint violation:" << endl;
            cout << ce.what() << endl;
          } catch (RuntimeException &re) {
            cout << "Runtime exception while generating CATEGORY level name for lipid " << lipid_name << endl;
          }
          slout << category << '\t';
          try {
            species = lipidAdduct->get_lipid_string(SPECIES);
          } catch (ConstraintViolationException &ce) {
            cout << "Constraint violation:" << endl;
            cout << ce.what() << endl;
          } catch (RuntimeException &re) {
            cout << "Runtime exception while generating SPECIES level name for lipid " << lipid_name << endl;
          }
          slout << species << endl;
          delete lipidAdduct;
          
          /* adding species into the map */
          if (category!="") {
            if (lipid_counts.find(category) == lipid_counts.end()) lipid_counts.insert({category, 0});
            ++lipid_counts.at(category);
          }
        } else {
          cout << "Could not parse '" << lipid_name << "'" << endl;
          slout << "N.A." << '\t' << "N.A." << '\t' << "N.A." << '\t' << "N.A." << endl;
        }
      }
      slout.close();
    }
    
    cout << "Parsed " << parsedLipids << " of " << totalLipids << " lipid names" << endl;    
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
