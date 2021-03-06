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


#ifndef FATTY_ACID_H
#define FATTY_ACID_H


#include <string>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidEnums.h"
#include <sstream>

using namespace std;
using namespace goslin;

class FattyAcid {
public:
    string name;
    int position;
    int num_carbon;
    int num_hydroxyl;
    int num_double_bonds;
    LipidFaBondType lipid_FA_bond_type;
    map<int, string> double_bond_positions;
    bool lcb;
    

    FattyAcid(string name, int num_carbon, int num_double_bonds, int num_hydroxyl, LipidFaBondType lipid_FA_bond_type, bool lcb, int position, map<int, string> *_double_bond_positions);
    FattyAcid();
    FattyAcid(FattyAcid* fa);
    string to_string(bool special_case = false);
    ElementTable* get_elements();
    static string suffix(LipidFaBondType _lipid_FA_bond_type);
    static bool lipid_FA_bond_type_prefix(LipidFaBondType lipid_FA_bond_type);
};
#endif /* FATTY_ACID_H */
