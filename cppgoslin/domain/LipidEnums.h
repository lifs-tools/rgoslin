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



#ifndef LIPID_ENUMS_H
#define LIPID_ENUMS_H


#include <vector>
#include <map>
#include <set>
#include "cppgoslin/domain/Element.h"


namespace goslin {
using namespace std;


enum LipidCategory {NO_CATEGORY,
    UNDEFINED,
    GL, // SLM:000117142 Glycerolipids
    GP, // SLM:000001193 Glycerophospholipids
    SP, // SLM:000000525 Sphingolipids
    ST, // SLM:000500463 Steroids and derivatives
    FA, // SLM:000390054 Fatty acyls and derivatives
    PK, // polyketides
    SL // Saccharo lipids
};


static const map<LipidCategory, string> CategoryString = {
    {NO_CATEGORY, "NO_CATEGORY"},
    {UNDEFINED, "UNDEFINED"},
    {GL, "GL"},
    {GP, "GP"},
    {SP, "SP"},
    {ST, "ST"},
    {FA, "FA"},
    {SL, "SL"}
};

    
#define is_level(x, y) (((x) & (y)) != 0)


enum LipidLevel {NO_LEVEL = 1,
    UNDEFINED_LEVEL = 2,
    CATEGORY = 4, // Mediators, Glycerolipids, Glycerophospholipids, Sphingolipids, Steroids, Prenols
    CLASS = 8, // Glyerophospholipids -> Glycerophosphoinositols PI
    SPECIES = 16, // Phosphatidylinositol 16:0 or PI 16:2;O2
    MOLECULAR_SPECIES = 32, // Phosphatidylinositol 8:2_8:0 or PI 8:2_8:0;O2
    SN_POSITION = 64, // Phosphatidylinositol 8:2/8:0 or PI 8:2/8:0;O2
    STRUCTURE_DEFINED = 128, // Phosphatidylinositol 8:2(4,6)/8:0 or PI 8:2(4,6)/8:0;(OH)2
    FULL_STRUCTURE = 256, // e.g. Phosphatidylinositol PI 8:2(4Z,6E)/8:2;3OH,7OH
    COMPLETE_STRUCTURE = 512, // e.g. Phosphatidylinositol PI 8:2(4Z,6E)[S]/8:2;3OH,7OH[R]
};



struct LipidClassMeta {
    LipidCategory lipid_category;
    string class_name;
    string description;
    int max_num_fa;
    int possible_num_fa;
    set<string> special_cases;
    ElementTable elements;
    vector<string> synonyms;
};


#include "cppgoslin/domain/ClassesEnum.h"


typedef map<LipidClass, LipidClassMeta> ClassMap;
enum LipidFaBondType {LCB_REGULAR = 0, LCB_EXCEPTION = 1, ETHER_PLASMANYL = 2, ETHER_PLASMENYL = 3, ETHER = 4, ETHER_UNSPECIFIED = 5, ESTER = 6, AMIDE = 7, UNDEFINED_FA = 8, NO_FA = 9};


static const set<LipidFaBondType> LCB_STATES {LCB_REGULAR, LCB_EXCEPTION};

class LipidClasses {
    public:
        static LipidClasses& get_instance()
        {
            static LipidClasses instance;
            return instance;
        }
    private:
        LipidClasses();
        
    public:
        ClassMap lipid_classes;
        LipidClasses(LipidClasses const&) = delete;
        void operator=(LipidClasses const&) = delete;
};



static map<LipidClass, string> ClassString;
static map<string, LipidClass> StringClass;
static map<string, LipidCategory> StringCategory;
}

#endif /* LIPID_ENUMS_H */
