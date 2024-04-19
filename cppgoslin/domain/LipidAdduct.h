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


#ifndef LIPID_ADDUCT_H
#define LIPID_ADDUCT_H

#include "cppgoslin/domain/FattyAcid.h"
#include <string>
#include <math.h>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/Element.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/LipidFullStructure.h"
#include "cppgoslin/domain/LipidStructureDefined.h"
#include "cppgoslin/domain/LipidSnPosition.h"
#include "cppgoslin/domain/LipidMolecularSpecies.h"
#include "cppgoslin/domain/LipidSpecies.h"
#include "cppgoslin/domain/Adduct.h"
#include <sstream>

using namespace std;
using namespace goslin;

class LipidAdduct {
public:
    LipidSpecies *lipid;
    Adduct *adduct;
    string sum_formula;
    
    LipidAdduct();
    LipidAdduct(LipidAdduct *lipid_adduct);
    ~LipidAdduct();
    string get_lipid_string(LipidLevel level = NO_LEVEL);
    string get_lipid_fragment_string(LipidLevel level = NO_LEVEL);
    string get_class_name();
    string get_extended_class();
    double get_mass();
    string get_sum_formula();
    LipidLevel get_lipid_level();
    ElementTable* get_elements();
    void sort_fatty_acyl_chains();
};

#endif /* LIPID_ADDUCT_H */
