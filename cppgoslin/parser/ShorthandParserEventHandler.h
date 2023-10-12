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


#ifndef SHORTHAND_PARSER_EVENT_HANDLER_H
#define SHORTHAND_PARSER_EVENT_HANDLER_H

#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/Adduct.h"
#include "cppgoslin/domain/Cycle.h"
#include "cppgoslin/domain/LipidAdduct.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/Headgroup.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/domain/GenericDatastructures.h"
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

class ShorthandParserEventHandler : public LipidBaseParserEventHandler {
public:
    vector<FunctionalGroup*> current_fas;
    GenericDictionary tmp;
    static const set<string> special_types;
    bool acer_species;
    bool contains_stereo_information;
    Element heavy_element;
    int heavy_element_number;
        
    ShorthandParserEventHandler();
    ~ShorthandParserEventHandler();
    void reset_lipid(TreeNode *node);
    void build_lipid(TreeNode *node);
    void add_cycle_element(TreeNode *node);
    void set_headgroup_name(TreeNode *node);
    void set_carbohydrate(TreeNode *node);
    void set_carbohydrate_sn_position(TreeNode *node);
    void set_carbohydrate_isomeric(TreeNode *node);
    void suffix_decorator_molecular(TreeNode *node);
    void suffix_decorator_species(TreeNode *node);
    void set_pl_hg_triple(TreeNode *node);
    void pre_sphingolipid(TreeNode *node);
    void set_ring_stereo(TreeNode *node);
    void post_sphingolipid(TreeNode *node);
    void set_hydroxyl(TreeNode *node);
    void new_lcb(TreeNode *node);
    void set_fatty_acyl_stereo(TreeNode *node);
    void add_pl_species_data(TreeNode *node);
    void new_fatty_acyl_chain(TreeNode *node);
    void add_fatty_acyl_chain(TreeNode *node);
    void add_dihydroxyl(TreeNode *node);
    void set_double_bond_position(TreeNode *node);
    void set_double_bond_information(TreeNode *node);
    void add_double_bond_information(TreeNode *node);
    void set_cistrans(TreeNode *node);
    void set_functional_group(TreeNode *node);
    void set_cycle(TreeNode *node);
    void add_cycle(TreeNode *node);
    void set_fatty_linkage_number(TreeNode *node);
    void set_hg_acyl(TreeNode *node);
    void add_hg_acyl(TreeNode *node);
    void set_hg_alkyl(TreeNode *node);
    void add_hg_alkyl(TreeNode *node);
    void set_linkage_type(TreeNode *node);
    void set_hydrocarbon_chain(TreeNode *node);
    void add_hydrocarbon_chain(TreeNode *node);
    void set_acyl_linkage(TreeNode *node);
    void add_acyl_linkage(TreeNode *node);
    void set_alkyl_linkage(TreeNode *node);
    void add_alkyl_linkage(TreeNode *node);
    void set_cycle_start(TreeNode *node);
    void set_cycle_end(TreeNode *node);
    void set_cycle_number(TreeNode *node);
    void set_cycle_db_count(TreeNode *node);
    void set_cycle_db_positions(TreeNode *node);
    void check_cycle_db_positions(TreeNode *node);
    void set_cycle_db_position(TreeNode *node);
    void set_cycle_db_position_cistrans(TreeNode *node);
    void set_functional_group_position(TreeNode *node);
    void set_functional_group_name(TreeNode *node);
    void set_functional_group_count(TreeNode *node);
    void set_functional_group_stereo(TreeNode *node);
    void set_sn_position_func_group(TreeNode *node);
    void add_functional_group(TreeNode *node);
    void set_ether_type(TreeNode *node);
    void set_ether_num(TreeNode *node);
    void set_species_level(TreeNode *node);
    void set_molecular_level(TreeNode *node);
    void set_carbon(TreeNode *node);
    void set_double_bond_count(TreeNode *node);
    void new_adduct(TreeNode *node);
    void add_adduct(TreeNode *node);
    void add_charge(TreeNode *node);
    void add_charge_sign(TreeNode *node);
    void set_acer(TreeNode *node);
    void set_acer_species(TreeNode *node);
    void set_sterol_definition(TreeNode *node);
    void set_carbohydrate_number(TreeNode *node);
    void set_glyco_sphingo_lipid(TreeNode *node);
    void set_heavy_element(TreeNode *node);
    void set_heavy_number(TreeNode *node);
    void add_heavy_component(TreeNode *node);
};



#endif /* SHORTHAND_PARSER_EVENT_HANDLER_H */
        
