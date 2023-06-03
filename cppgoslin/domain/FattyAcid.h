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


#ifndef FATTY_ACID_H
#define FATTY_ACID_H


#include <string>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/domain/DoubleBonds.h"
#include "cppgoslin/domain/StringFunctions.h"
#include <sstream>

using namespace std;
using namespace goslin;

class FattyAcid : public FunctionalGroup {
public:
    int num_carbon;
    LipidFaBondType lipid_FA_bond_type;
    
    FattyAcid(string name, int num_carbon = 0, DoubleBonds* double_bonds = 0, map<string, vector<FunctionalGroup*> >* functional_groups = 0, LipidFaBondType lipid_FA_bond_type = ESTER, int position = 0);
    string to_string(LipidLevel level);
    void compute_elements();
    FattyAcid* copy();
    int get_double_bonds();
    void set_type(LipidFaBondType _lipid_FA_bond_type);
    ElementTable* get_functional_group_elements();
    static string get_prefix(LipidFaBondType _lipid_FA_bond_type);
    int num_oxygens();
    static bool lipid_FA_bond_type_prefix(LipidFaBondType lipid_FA_bond_type);
    const set<string> fg_exceptions {"acyl", "alkyl", "cy", "cc", "acetoxy"};
};


    


class AcylAlkylGroup : public FunctionalGroup {
public:
    bool alkyl;
    bool N_bond;
    
    
    AcylAlkylGroup(FattyAcid* _fa, int _position = -1, int count = 1, bool _alkyl = false, bool N_bond = false);
    AcylAlkylGroup* copy();
    void set_N_bond_type(bool N_bond);
    string to_string(LipidLevel level);
};
    

    
class CarbonChain : public FunctionalGroup {
public:
    CarbonChain(FattyAcid* _fa, int _position = -1, int _count = 1);
    CarbonChain* copy();
    string to_string(LipidLevel level);
};

#endif /* FATTY_ACID_H */
