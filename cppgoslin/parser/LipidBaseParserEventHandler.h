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


#ifndef LIPID_BASE_PARSER_EVENT_HANDLER_H
#define LIPID_BASE_PARSER_EVENT_HANDLER_H

#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/Adduct.h"
#include "cppgoslin/domain/LipidAdduct.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/Cycle.h"
#include "cppgoslin/domain/Headgroup.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/parser/BaseParserEventHandler.h"
#include "cppgoslin/parser/Parser.h"
#include <string>
#include <math.h>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
using namespace goslin;

class LipidBaseParserEventHandler : public BaseParserEventHandler<LipidAdduct*> {
public:
    LipidLevel level;
    string head_group;
    FattyAcid *lcb;
    vector<FattyAcid*> *fa_list;
    FattyAcid *current_fa;
    vector<HeadgroupDecorator*> *headgroup_decorators;
    bool use_head_group;
    static const set<string> SP_EXCEPTION_CLASSES;
    Adduct* adduct;
    static const map<string, int> fa_synonyms;
        
    LipidBaseParserEventHandler();
    ~LipidBaseParserEventHandler();
    void set_lipid_level(LipidLevel _level);
    bool sp_regular_lcb();
    Headgroup* prepare_headgroup_and_checks(bool allow_class_shift = true);
    LipidSpecies *assemble_lipid(Headgroup *headgroup);
    FattyAcid* resolve_fa_synonym(string mediator_name);
    bool check_full_structure(FunctionalGroup *obj);
};


#endif /* LIPID_BASE_PARSER_EVENT_HANDLER_H */
        
