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


#ifndef LIPID_MAPS_PARSER_EVENT_HANDLER_H
#define LIPID_MAPS_PARSER_EVENT_HANDLER_H

#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/Adduct.h"
#include "cppgoslin/domain/LipidAdduct.h"
#include "cppgoslin/domain/Cycle.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/domain/Headgroup.h"
#include "cppgoslin/parser/LipidBaseParserEventHandler.h"
#include <string>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
using namespace goslin;

static const set<string> head_group_exceptions = {"PA", "PC", "PE", "PG", "PI", "PS"};

class LipidMapsParserEventHandler : public LipidBaseParserEventHandler {
public:
    bool omit_fa;
    int db_numbers;
    int db_position;
    string db_cistrans;
    string mod_text;
    int mod_pos;
    int mod_num;
    bool add_omega_linoleoyloxy_Cer;
    

    LipidMapsParserEventHandler();
    ~LipidMapsParserEventHandler();
    void reset_lipid(TreeNode* node);
    void set_molecular_subspecies_level(TreeNode* node);
    void mediator_event(TreeNode* node);
    void set_head_group_name(TreeNode* node);
    void set_species_level(TreeNode* node);
    void add_functional_group(TreeNode* node);
    void set_structural_subspecies_level(TreeNode* node);
    void set_mod(TreeNode* node);
    void set_mod_text(TreeNode* node);
    void set_mod_pos(TreeNode* node);
    void set_mod_num(TreeNode* node);
    void new_fa(TreeNode* node);
    void new_lcb(TreeNode* node);
    void clean_lcb(TreeNode* node);
    void append_fa(TreeNode* node);
    void add_ether(TreeNode* node);
    void set_sphingoxine(TreeNode* node);
    void add_hydroxyl(TreeNode* node);
    void add_double_bonds(TreeNode* node);
    void add_carbon(TreeNode* node);
    void build_lipid(TreeNode* node);
    void add_hydroxyl_lcb(TreeNode* node);
    void pure_fa(TreeNode* node);
    void add_glyco(TreeNode* node);
    void set_isomeric_level(TreeNode* node);
    void add_db_position(TreeNode* node);
    void add_db_position_number(TreeNode* node);
    void add_cistrans(TreeNode* node);
    void set_omega_head_group_name(TreeNode* node);
    void add_ACer(TreeNode* node);
    void new_adduct(TreeNode *node);
    void add_adduct(TreeNode *node);
    void add_charge(TreeNode *node);
    void add_charge_sign(TreeNode *node);
    
    static const map<string, int> acer_heads;
        
};



#endif /* LIPID_MAPS_PARSER_EVENT_HANDLER_H */
