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


#ifndef GOSLIN_PARSER_EVENT_HANDLER_H
#define GOSLIN_PARSER_EVENT_HANDLER_H

#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/Adduct.h"
#include "cppgoslin/domain/LipidAdduct.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/Headgroup.h"
#include "cppgoslin/domain/Cycle.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/parser/LipidBaseParserEventHandler.h"
#include <string>
#include <math.h>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
using namespace goslin;

class GoslinParserEventHandler : public LipidBaseParserEventHandler {
public:
    int db_position;
    string db_cistrans;
    bool unspecified_ether;
    char plasmalogen;
    string mediator_function;
    vector<int> mediator_function_positions;
    bool mediator_suffix;
    Element heavy_element;
    int heavy_element_number;
    
    static const map<string, int> mediator_FA;
    static const map<string, int> mediator_DB;
        
    GoslinParserEventHandler();
    ~GoslinParserEventHandler();
    void reset_lipid(TreeNode *node);
    void set_head_group_name(TreeNode *node);
    void set_species_level(TreeNode *node);
    void set_molecular_subspecies_level(TreeNode *node);
    void new_fa(TreeNode *node);
    void new_lcb(TreeNode *node);
    void clean_lcb(TreeNode *node);
    void append_fa(TreeNode *node);
    void build_lipid(TreeNode *node);
    void add_ether(TreeNode *node);
    void add_old_hydroxyl(TreeNode *node);
    void add_double_bonds(TreeNode *node);
    void add_carbon(TreeNode *node);
    void add_hydroxyl(TreeNode *node);
    void new_adduct(TreeNode *node);
    void add_adduct(TreeNode *node);
    void add_charge(TreeNode *node);
    void add_charge_sign(TreeNode *node);
    void set_unspecified_ether(TreeNode *node);
    void set_plasmalogen(TreeNode *node);
    
    void set_isomeric_level(TreeNode* node);
    void add_db_position(TreeNode* node);
    void add_db_position_number(TreeNode* node);
    void add_cistrans(TreeNode* node);
    
    void set_mediator(TreeNode *node);
    void set_unstructured_mediator(TreeNode *node);
    void set_trivial_mediator(TreeNode *node);
    void set_mediator_carbon(TreeNode *node);
    void set_mediator_db(TreeNode *node);
    void set_mediator_function(TreeNode *node);
    void set_mediator_function_position(TreeNode *node);
    void add_mediator_function(TreeNode *node);
    void add_mediator_suffix(TreeNode *node);
    void add_mediator(TreeNode *node);
    void set_mediator_tetranor(TreeNode *node);
    
    void set_heavy_d_element(TreeNode *node);
    void set_heavy_d_number(TreeNode *node);
    void set_heavy_element(TreeNode *node);
    void set_heavy_number(TreeNode *node);
    void add_heavy_component(TreeNode *node);
};


#endif /* GOSLIN_PARSER_EVENT_HANDLER_H */
        
