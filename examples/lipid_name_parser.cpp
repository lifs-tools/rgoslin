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


#include "cppgoslin/cppgoslin.h"
#include <iostream>
#include <string>

using namespace std;
using namespace goslin;

int main(){
    /* create instance of lipid parser containing several grammars */
    LipidParser lipid_parser;
    
    /* parsing lipid name into a lipid container data structure */
    string lipid_name = "PA(12:0(2OH)/14:1(3Z))";
    LipidAdduct* lipid = lipid_parser.parse(lipid_name);
    
    /* checking if parsing was successful, otherwise lipid reference remains NULL */
    if (lipid != NULL){
        /* creating lipid name according to the recent lipid nomenclature */
        cout << "input:\t\t\t\t" << lipid_name << endl;
        cout << "Full structure level:\t\t" << lipid->get_lipid_string(FULL_STRUCTURE) << endl;
        cout << "Structure defined level:\t" << lipid->get_lipid_string(STRUCTURE_DEFINED) << endl;
        cout << "sn-position level:\t\t" << lipid->get_lipid_string(SN_POSITION) << endl;
        cout << "Molecular species level:\t" << lipid->get_lipid_string(MOLECULAR_SPECIES) << endl;
        cout << "Species level:\t\t\t" << lipid->get_lipid_string(SPECIES) << endl;
        cout << "Class name:\t\t\t" << lipid->get_class_name() << endl;
        
        /* important to delete lipid to avoid memory leaks */
        delete lipid;
    }
    
    return 0;
}
