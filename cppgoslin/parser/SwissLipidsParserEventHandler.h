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


#ifndef SWISS_LIPIDS_PARSER_EVENT_HANDLER_H
#define SWISS_LIPIDS_PARSER_EVENT_HANDLER_H

#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/Adduct.h"
#include "cppgoslin/domain/StringFunctions.h"
#include "cppgoslin/domain/LipidAdduct.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/Headgroup.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/parser/LipidBaseParserEventHandler.h"
#include <string>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
using namespace goslin;

class SwissLipidsParserEventHandler : public LipidBaseParserEventHandler {
public:
    int db_position;
    string db_cistrans;
    int suffix_number;
        
    SwissLipidsParserEventHandler();
    ~SwissLipidsParserEventHandler();
    
    void reset_lipid(TreeNode *node);
    void build_lipid(TreeNode *node);
    void set_head_group_name(TreeNode *node);
    void set_species_level(TreeNode *node);
    void set_molecular_level(TreeNode *node);
    void set_level(LipidLevel _level);
    void new_lcb(TreeNode *node);
    void clean_lcb(TreeNode *node);
    void new_fa(TreeNode *node);
    void append_fa(TreeNode *node);
    void add_ether(TreeNode *node);
    void add_hydroxyl(TreeNode *node);
    void add_double_bonds(TreeNode *node);
    void add_carbon(TreeNode *node);
    void mediator_event(TreeNode* node);
    void set_nape(TreeNode *node);
    void set_isomeric_level(TreeNode* node);
    void add_db_position(TreeNode* node);
    void add_db_position_number(TreeNode* node);
    void add_cistrans(TreeNode* node);
    void set_species_fa(TreeNode *node);
    void set_head_group_name_se(TreeNode *node);
    void add_one_hydroxyl(TreeNode *node);
    void add_suffix_number(TreeNode *node);
    void add_fa_lcb_suffix_type(TreeNode *node);
    void new_adduct(TreeNode *node);
    void add_adduct(TreeNode *node);
    void add_charge(TreeNode *node);
    void add_charge_sign(TreeNode *node);
    
};


#endif /* SWISS_LIPIDS_PARSER_EVENT_HANDLER_H */
        
