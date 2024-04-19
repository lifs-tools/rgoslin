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


#include "cppgoslin/parser/SwissLipidsParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<LipidAdduct*>::registered_events->insert({x, bind(&SwissLipidsParserEventHandler::y, this, placeholders::_1)})
    

SwissLipidsParserEventHandler::SwissLipidsParserEventHandler() : LipidBaseParserEventHandler() {
    
    reg("lipid_pre_event", reset_lipid);
    reg("lipid_post_event", build_lipid);
    reg("fa_hg_pre_event", set_head_group_name);
    reg("gl_hg_pre_event", set_head_group_name);
    reg("gl_molecular_hg_pre_event", set_head_group_name);
    reg("mediator_pre_event", mediator_event);
    reg("gl_mono_hg_pre_event", set_head_group_name);
    reg("pl_hg_pre_event", set_head_group_name);
    reg("pl_three_hg_pre_event", set_head_group_name);
    reg("pl_four_hg_pre_event", set_head_group_name);
    reg("sl_hg_pre_event", set_head_group_name);
    reg("st_species_hg_pre_event", set_head_group_name);
    reg("st_sub1_hg_pre_event", set_head_group_name);
    reg("st_sub2_hg_pre_event", set_head_group_name_se);
    reg("fa_species_pre_event", set_species_level);
    reg("gl_molecular_pre_event", set_molecular_level);
    reg("unsorted_fa_separator_pre_event", set_molecular_level);
    reg("fa2_unsorted_pre_event", set_molecular_level);
    reg("fa3_unsorted_pre_event", set_molecular_level);
    reg("fa4_unsorted_pre_event", set_molecular_level);
    reg("db_single_position_pre_event", set_isomeric_level);
    reg("db_single_position_post_event", add_db_position);
    reg("db_position_number_pre_event", add_db_position_number);
    reg("cistrans_pre_event", add_cistrans);
    reg("lcb_pre_event", new_lcb);
    reg("lcb_post_event", clean_lcb);
    reg("fa_pre_event", new_fa);
    reg("fa_post_event", append_fa);
    reg("ether_pre_event", add_ether);
    reg("hydroxyl_pre_event", add_hydroxyl);
    reg("db_count_pre_event", add_double_bonds);
    reg("carbon_pre_event", add_carbon);
    reg("sl_lcb_species_pre_event", set_species_level);
    reg("st_species_fa_post_event", set_species_fa);
    reg("fa_lcb_suffix_type_pre_event", add_fa_lcb_suffix_type);
    reg("fa_lcb_suffix_number_pre_event", add_suffix_number);
    reg("pl_three_post_event", set_nape);
    reg("adduct_info_pre_event", new_adduct);
    reg("adduct_pre_event", add_adduct);
    reg("charge_pre_event", add_charge);
    reg("charge_sign_pre_event", add_charge_sign);
    
    debug = "";
}


SwissLipidsParserEventHandler::~SwissLipidsParserEventHandler(){
}


void SwissLipidsParserEventHandler::reset_lipid(TreeNode *node) {
    level = FULL_STRUCTURE;
    head_group = "";
    lcb = NULL;
    fa_list->clear();
    current_fa = NULL;
    adduct = 0;
    use_head_group = false;
    db_position = 0;
    db_cistrans = "";
    headgroup_decorators->clear();
    suffix_number = -1;
}


void SwissLipidsParserEventHandler::set_nape(TreeNode *node){
    head_group = "PE-N";
    HeadgroupDecorator* hgd = new HeadgroupDecorator("decorator_acyl", -1, 1, 0, true);
    headgroup_decorators->push_back(hgd);
    hgd->functional_groups->insert({"decorator_acyl", vector<FunctionalGroup*>()});
    hgd->functional_groups->at("decorator_acyl").push_back(fa_list->at(fa_list->size() - 1));
    fa_list->pop_back();
}


void SwissLipidsParserEventHandler::set_isomeric_level(TreeNode* node){
    db_position = 0;
    db_cistrans = "";
}


void SwissLipidsParserEventHandler::add_db_position(TreeNode* node){
    if (current_fa != NULL){
        current_fa->double_bonds->double_bond_positions.insert({db_position, db_cistrans});
        if (db_cistrans != "E" && db_cistrans != "Z") set_lipid_level(STRUCTURE_DEFINED);
    }
}


void SwissLipidsParserEventHandler::add_db_position_number(TreeNode* node){
    db_position = node->get_int();
}


void SwissLipidsParserEventHandler::add_cistrans(TreeNode* node){
    db_cistrans = node->get_text();
}


void SwissLipidsParserEventHandler::set_head_group_name(TreeNode *node) {
    head_group = node->get_text();
}


void SwissLipidsParserEventHandler::set_head_group_name_se(TreeNode *node){
    head_group = replace_all(node->get_text(), "(", " ");
}



void SwissLipidsParserEventHandler::set_species_level(TreeNode *node) {
    set_lipid_level(SPECIES);
}
    



void SwissLipidsParserEventHandler::set_molecular_level(TreeNode *node) {
    set_lipid_level(MOLECULAR_SPECIES);
}


void SwissLipidsParserEventHandler::mediator_event(TreeNode* node){
    use_head_group = true;
    head_group = node->get_text();
}
    
    

void SwissLipidsParserEventHandler::new_fa(TreeNode *node) {
    current_fa = new FattyAcid("FA");
}
    
    

void SwissLipidsParserEventHandler::new_lcb(TreeNode *node) {
    lcb = new FattyAcid("LCB");
    lcb->set_type(LCB_REGULAR);
    current_fa = lcb;
    set_lipid_level(STRUCTURE_DEFINED);
}
        
        

void SwissLipidsParserEventHandler::clean_lcb(TreeNode *node) {
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
    current_fa = NULL;
}
    
    
        

void SwissLipidsParserEventHandler::append_fa(TreeNode *node) {
    if (current_fa->double_bonds->get_num() < 0){
        throw LipidException("Double bond count does not match with number of double bond positions");
    }
    
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
    
    fa_list->push_back(current_fa);
    current_fa = NULL;
}
    
    

void SwissLipidsParserEventHandler::build_lipid(TreeNode *node) {
    if (lcb){
        set_lipid_level(STRUCTURE_DEFINED);
        fa_list->insert(fa_list->begin(), lcb);
    }

    Headgroup *headgroup = prepare_headgroup_and_checks();
    
    LipidAdduct *lipid = new LipidAdduct();
    lipid->lipid = assemble_lipid(headgroup);
    lipid->adduct = adduct;
    BaseParserEventHandler<LipidAdduct*>::content = lipid;
}
    
    

void SwissLipidsParserEventHandler::add_ether(TreeNode *node) {
    string ether = node->get_text();
    if (ether == "O-") current_fa->lipid_FA_bond_type = ETHER_PLASMANYL;
    else if (ether == "P-") current_fa->lipid_FA_bond_type = ETHER_PLASMENYL;
}
    
    

void SwissLipidsParserEventHandler::add_hydroxyl(TreeNode *node) {
    string old_hydroxyl = node->get_text();
    int num_h = 0;
    if (old_hydroxyl == "m") num_h = 1;
    else if (old_hydroxyl == "d") num_h = 2;
    else if (old_hydroxyl == "t") num_h = 3;
    
    if (sp_regular_lcb()) num_h -= 1;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("OH").push_back(functional_group);
}


void SwissLipidsParserEventHandler::add_one_hydroxyl(TreeNode *node) {
    if (contains_val_p(current_fa->functional_groups, "OH") && current_fa->functional_groups->at("OH").at(0)->position == -1){
        current_fa->functional_groups->at("OH").at(0)->count += 1;
    }
    else {
        FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
        if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
        current_fa->functional_groups->at("OH").push_back(functional_group);
    }
}
    
    

void SwissLipidsParserEventHandler::add_double_bonds(TreeNode *node) {
    current_fa->double_bonds->num_double_bonds += node->get_int();
}



void SwissLipidsParserEventHandler::add_suffix_number(TreeNode *node){
    suffix_number = node->get_int();
}



void SwissLipidsParserEventHandler::add_fa_lcb_suffix_type(TreeNode *node){
    string suffix_type = node->get_text();
    if (suffix_type == "me"){
        suffix_type = "Me";
        current_fa->num_carbon -= 1;
    }
        
    FunctionalGroup *functional_group = KnownFunctionalGroups::get_functional_group(suffix_type);
    functional_group->position = suffix_number;
    if (functional_group->position == -1) set_lipid_level(STRUCTURE_DEFINED);
    if (uncontains_val_p(current_fa->functional_groups, suffix_type)) current_fa->functional_groups->insert({suffix_type, vector<FunctionalGroup*>()});
    current_fa->functional_groups->at(suffix_type).push_back(functional_group);
            
    suffix_number = -1;
}
    
    

void SwissLipidsParserEventHandler::add_carbon(TreeNode *node) {
    current_fa->num_carbon = node->get_int();
}

        

void SwissLipidsParserEventHandler::set_species_fa(TreeNode *node){
    head_group += " 27:1";
    fa_list->at(fa_list->size() -1)->num_carbon -= 27;
    fa_list->at(fa_list->size() -1)->double_bonds->num_double_bonds -= 1;
}
    
    

void SwissLipidsParserEventHandler::new_adduct(TreeNode *node) {
    if (!adduct) adduct = new Adduct("", "");
}
    
    

void SwissLipidsParserEventHandler::add_adduct(TreeNode *node) {
    adduct->adduct_string = node->get_text();
}
    
    

void SwissLipidsParserEventHandler::add_charge(TreeNode *node) {
    adduct->charge = node->get_int();
}
    
    

void SwissLipidsParserEventHandler::add_charge_sign(TreeNode *node) {
    string sign = node->get_text();
    if (sign == "+") adduct->set_charge_sign(1);
    else if (sign == "-") adduct->set_charge_sign(-1);
    if (adduct->charge == 0) adduct->charge = 1;
}
        
