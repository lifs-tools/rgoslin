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


#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>
#include "cppgoslin/domain/StringFunctions.h"
#include "cppgoslin/parser/Parser.h"
#include "cppgoslin/parser/SumFormulaParserEventHandler.h"

using namespace std;
using namespace goslin;



void writeLipidEnum(string ofFileName){
    ofstream offile(ofFileName.c_str());
    
    ifstream infile("data/goslin/lipid-list.csv");
    if (!infile.good()){
        cout << "Error: file 'data/goslin/lipid-list.csv' not found." << endl;
        exit(-1);
    }
    
    ifstream functional_file("data/goslin/functional-groups.csv");
    if (!functional_file.good()){
        cout << "Error: file 'data/goslin/lipid-list.csv' not found." << endl;
        exit(-1);
    }
    
    string line;
    unsigned int i = 0;
    int SYNONYM_START_INDEX = 7;
    map<string, int> enum_names = {{"GL", 1}, {"GP", 1}, {"SP", 1}, {"ST", 1}, {"FA", 1}, {"PK", 1}, {"SL", 1}, {"UNDEFINED", 1}};
    
    map<string, vector<string>*> data;
    vector< vector<string>*> functional_data;
    set<string> functional_data_set;
    set<string> keys;
    while (getline(infile, line)){
        if (i++ == 0) continue;
        vector<string>* tokens = split_string(line, ',', '"', true);
        for (unsigned int i = 0; i < tokens->size(); ++i){
            string s = tokens->at(i);
            if (s.length() >= 2 && s[0] == '"' && s[s.length() - 1] == '"'){
                tokens->at(i) = s.substr(1, s.length() - 2);
            }
        }
        
        if (keys.find(tokens->at(0)) != keys.end()){
            cout << "Error: lipid name '" << tokens->at(0) << "' occurs multiple times in the lipid list." << endl;
            exit(-1);
        }
        keys.insert(tokens->at(0));
        
        
        for (int i = SYNONYM_START_INDEX; i < (int)tokens->size(); ++i){
            string test_lipid_name = tokens->at(i);
            if (test_lipid_name == "") continue;
            if (keys.find(test_lipid_name) != keys.end()){
                cout << "Error: lipid name '" << test_lipid_name << "' occurs multiple times in the lipid list." << endl;
                
                exit(-1);
            }
            keys.insert(test_lipid_name);
        }
        
        
        string enum_name = tokens->at(0);
        
        for (unsigned int i = 0; i < enum_name.length(); ++i){
            char c = enum_name[i];
            if ('A' <= c && c <= 'Z'){
                
            }
            else if ('0' <= c && c <= '9'){
                
            }
            else if ('a' <= c && c <= 'z'){
                enum_name[i] = c - ('a' - 'A');
            }
            else {
                enum_name[i] = '_';
            }
        }
        
        
        if (enum_name[0] == '_'){
            enum_name = "L" + enum_name;
        }
        
        if (enum_name[0] < 'A' || 'Z' < enum_name[0]){
            enum_name = "L" + enum_name;
        }
        
        if (enum_names.find(enum_name) == enum_names.end()){
            enum_names.insert({enum_name, 1});
        }
        else {
            int cnt = enum_names.at(enum_name)++;
            enum_names.at(enum_name) = cnt;
            enum_name += ('A' + cnt - 1);
            enum_names.insert({enum_name, 1});
        }
        
        data.insert({enum_name, tokens});
    }
    SumFormulaParserEventHandler sum_formula_handler;
    Parser<ElementTable*> parser(&sum_formula_handler, "data/goslin/SumFormula.g4", DEFAULT_QUOTE);
    
    while (getline(functional_file, line)){
        if (i++ == 0) continue;
        functional_data.push_back(split_string(line, ',', '"', true));
        string fd_name = functional_data.back()->at(1);
        if (functional_data_set.find(fd_name) != functional_data_set.end()){
            cout << "Error: functional group '" + fd_name + "' occurs multiple times!" << endl;
            exit(-1);
        }
        functional_data_set.insert(fd_name);
    }
        
    
    
    ofstream enums("cppgoslin/domain/ClassesEnum.h");
    enums << "enum LipidClass {NO_CLASS, UNDEFINED_CLASS";
    
    
     map<Element, string> table_symbol{{ELEMENT_C, "ELEMENT_C"}, {ELEMENT_H, "ELEMENT_H"}, {ELEMENT_N, "ELEMENT_N"}, {ELEMENT_O, "ELEMENT_O"}, {ELEMENT_P, "ELEMENT_P"}, {ELEMENT_S, "ELEMENT_S"}, {ELEMENT_H2, "ELEMENT_H2"}, {ELEMENT_C13, "ELEMENT_C13"}, {ELEMENT_N15, "ELEMENT_N15"}, {ELEMENT_O17, "ELEMENT_O17"}, {ELEMENT_O18, "ELEMENT_O18"}, {ELEMENT_P32, "ELEMENT_P32"}, {ELEMENT_S33, "ELEMENT_S33"}, {ELEMENT_S34, "ELEMENT_S34"}, {ELEMENT_F, "ELEMENT_F"}, {ELEMENT_Cl, "ELEMENT_Cl"}, {ELEMENT_Br, "ELEMENT_Br"}, {ELEMENT_I, "ELEMENT_I"}, {ELEMENT_As, "ELEMENT_As"}};
     
    
    
    offile << "/* DO NOT CHANGE THE FILE, IT IS AUTOMATICALLY GENERATED */" << endl << endl;
    
    
    offile << "/*" << endl;
    offile << "MIT License" << endl;
    offile << endl;
    offile << "Copyright (c) the authors (listed in global LICENSE file)" << endl;
    offile << endl;
    offile << "Permission is hereby granted, free of charge, to any person obtaining a copy" << endl;
    offile << "of this software and associated documentation files (the \"Software\"), to deal" << endl;
    offile << "in the Software without restriction, including without limitation the rights" << endl;
    offile << "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell" << endl;
    offile << "copies of the Software, and to permit persons to whom the Software is" << endl;
    offile << "furnished to do so, subject to the following conditions:" << endl;
    offile << endl;
    offile << "The above copyright notice and this permission notice shall be included in all" << endl;
    offile << "copies or substantial portions of the Software." << endl;
    offile << endl;
    offile << "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR" << endl;
    offile << "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY," << endl;
    offile << "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE" << endl;
    offile << "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER" << endl;
    offile << "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM," << endl;
    offile << "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE" << endl;
    offile << "SOFTWARE." << endl;
    offile << "*/" << endl;
    offile << endl;
    offile << "#include \"cppgoslin/domain/LipidEnums.h\"" << endl;
    offile << "#include \"cppgoslin/domain/FunctionalGroup.h\"" << endl;
    offile << endl;
    offile << "using namespace std;" << endl;
    offile << "using namespace goslin;" << endl;
    offile << endl;
    offile << endl;
    offile << "LipidClasses::LipidClasses(){" << endl;
    offile << "    lipid_classes = {" << endl;
    unsigned int cnt = 0;
    for (auto& kv : data){
        // add lipid category, description, max num fa, possible num fa
        offile << "    {" << kv.first << ", {" << kv.second->at(1) << ", \"" << kv.second->at(0) << "\", \"" << kv.second->at(2) << "\", ";
        offile << kv.second->at(3) << ", " << kv.second->at(4) << ", {";
        enums << ", " << kv.first;
        
        vector<string>* tokens = split_string(kv.second->at(5), ';', '"');
        for (unsigned int i = 0; i < tokens->size(); ++i){
            string tok = strip(tokens->at(i), ' ');
            if (i > 0) offile << ", ";
            offile << "\"" << tok << "\"";
        }
        delete tokens;
        offile << "}, {";
        
        
        
        // add element table
        ElementTable* table = kv.second->at(6).length() > 0 ? parser.parse(kv.second->at(6)) : create_empty_table();
        int ii = 0;
        for (auto& table_kv : *table){
            if (ii++ > 0) offile << ", ";
            offile << "{" << table_symbol.at(table_kv.first) << ", " << table_kv.second << "}";
        }
        delete table;
         
        // add synonyms
        offile << "}, {\"" << kv.second->at(0) << "\"";
        for (unsigned int i = SYNONYM_START_INDEX; i < kv.second->size(); ++i){
            string synonym = kv.second->at(i);
            if (synonym.length() < 1) continue;
            offile << ", \"" << synonym << "\"";
        }
        offile << "} } }" << (++cnt < data.size() ? ",\n" : "\n") << endl;
    }
    
    enums << "};" << endl;

    
    offile << "    };" << endl; 
    offile << "}" << endl; 
    offile << endl;
    offile << endl;
    offile << endl;
    offile << endl;
    
    
    offile << "KnownFunctionalGroups::KnownFunctionalGroups(){" << endl;
    
    cnt = 0;
    for (auto &row : functional_data){
        if (cnt++ == 0) continue;
        
        row->push_back(row->at(1));
        for (int i = 7; i < (int)row->size(); ++i){
            if (row->at(i) == "") continue;
            
            offile << "    // " << row->at(6) << endl;
            if (row->at(0) == "FG"){
                offile << "    known_functional_groups.insert({\"" << row->at(i) << "\", new FunctionalGroup(\"" << row->at(1) << "\", -1, 1, new DoubleBonds(" << row->at(3) << "), " << row->at(4) << ", \"\", " << row->at(5) << ", new ElementTable{";
                
                // add element table
                ElementTable* table = row->at(2).length() > 0 ? parser.parse(row->at(2)) : create_empty_table();
                int ii = 0;
                for (auto& table_kv : *table){
                    if (ii++ > 0) offile << ", ";
                    offile << "{" << table_symbol.at(table_kv.first) << ", " << table_kv.second << "}";
                }
                delete table;
                offile << "}, 0)}";
            }
            else {
                offile << "    known_functional_groups.insert({\"" << row->at(i) << "\", new HeadgroupDecorator(\"" << row->at(1) << "\", -1, 1, new ElementTable{";
                
                // add element table
                ElementTable* table = row->at(2).length() > 0 ? parser.parse(row->at(2)) : create_empty_table();
                int ii = 0;
                for (auto& table_kv : *table){
                    if (ii++ > 0) offile << ", ";
                    offile << "{" << table_symbol.at(table_kv.first) << ", " << table_kv.second << "}";
                }
                delete table;
                offile << "})}";
            }
            offile << ");" << endl;
            offile << "\n\n";
        }
        
    }
    
    offile << "}" << endl;
    offile << endl;
    
    
    for (auto &f : functional_data) delete f;
}








int main(int argc, char** argv){
    if (argc < 2){
        cout << "Error, specify grammar output filename." << endl;
        exit(-1);
    }

    writeLipidEnum(argv[1]);
    return 0;
}
