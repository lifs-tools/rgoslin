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


#include "cppgoslin/parser/LipidMapsParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<LipidAdduct*>::registered_events->insert({x, bind(&LipidMapsParserEventHandler::y, this, placeholders::_1)})
    

LipidMapsParserEventHandler::LipidMapsParserEventHandler() : LipidBaseParserEventHandler() {
    reg("lipid_pre_event", reset_lipid);
    reg("lipid_post_event", build_lipid);
    
    reg("mediator_pre_event", mediator_event);
    
    reg("sgl_species_pre_event", set_species_level);
    reg("species_fa_pre_event", set_species_level);
    reg("tgl_species_pre_event", set_species_level);
    reg("dpl_species_pre_event", set_species_level);
    reg("cl_species_pre_event", set_species_level);
    reg("dsl_species_pre_event", set_species_level);
    reg("fa2_unsorted_pre_event", set_molecular_subspecies_level);
    reg("fa3_unsorted_pre_event", set_molecular_subspecies_level);
    reg("fa4_unsorted_pre_event", set_molecular_subspecies_level);
    reg("hg_dg_pre_event", set_molecular_subspecies_level);
    reg("fa_lpl_molecular_pre_event", set_molecular_subspecies_level);
    reg("hg_lbpa_pre_event", set_molecular_subspecies_level);

    reg("fa_no_hg_pre_event", pure_fa);
    
    reg("hg_sgl_pre_event", set_head_group_name);
    reg("hg_gl_pre_event", set_head_group_name);
    reg("hg_cl_pre_event", set_head_group_name);
    reg("hg_dpl_pre_event", set_head_group_name);
    reg("hg_lpl_pre_event", set_head_group_name);
    reg("hg_threepl_pre_event", set_head_group_name);
    reg("hg_fourpl_pre_event", set_head_group_name);
    reg("hg_dsl_pre_event", set_head_group_name);
    reg("hg_cpa_pre_event", set_head_group_name);
    reg("ch_pre_event", set_head_group_name);
    reg("hg_che_pre_event", set_head_group_name);
    reg("mediator_const_pre_event", set_head_group_name);
    reg("pk_hg_pre_event", set_head_group_name);
    reg("hg_fa_pre_event", set_head_group_name);
    reg("hg_lsl_pre_event", set_head_group_name);
    reg("special_cer_pre_event", set_head_group_name);
    reg("special_cer_hg_pre_event", set_head_group_name);
    reg("omega_linoleoyloxy_Cer_pre_event", set_omega_head_group_name);
    
    reg("lcb_pre_event", new_lcb);
    reg("lcb_post_event", clean_lcb);
    reg("fa_pre_event", new_fa);
    reg("fa_post_event", append_fa);
    
    reg("glyco_struct_pre_event", add_glyco);
    
    reg("db_single_position_pre_event", set_isomeric_level);
    reg("db_single_position_post_event", add_db_position);
    reg("db_position_number_pre_event", add_db_position_number);
    reg("cistrans_pre_event", add_cistrans);
    
    reg("ether_pre_event", add_ether);
    reg("hydroxyl_pre_event", add_hydroxyl);
    reg("hydroxyl_lcb_pre_event", add_hydroxyl_lcb);
    reg("db_count_pre_event", add_double_bonds);
    reg("carbon_pre_event", add_carbon);
    
    reg("structural_mod_pre_event", set_structural_subspecies_level);
    reg("single_mod_pre_event", set_mod);
    reg("mod_text_pre_event", set_mod_text);
    reg("mod_pos_pre_event", set_mod_pos);
    reg("mod_num_pre_event", set_mod_num);
    reg("single_mod_post_event", add_functional_group);
    reg("special_cer_prefix_pre_event", add_ACer);
    
    
    reg("adduct_info_pre_event", new_adduct);
    reg("adduct_pre_event", add_adduct);
    reg("charge_pre_event", add_charge);
    reg("charge_sign_pre_event", add_charge_sign);
    
    debug = "";
} 


LipidMapsParserEventHandler::~LipidMapsParserEventHandler(){
}
    
    
void LipidMapsParserEventHandler::reset_lipid(TreeNode* node){
    level = FULL_STRUCTURE;
    head_group = "";
    lcb = NULL;
    fa_list->clear();
    current_fa = NULL;
    use_head_group = false;
    omit_fa = false;
    db_position = 0;
    adduct = 0;
    db_numbers = -1;
    db_cistrans = "";
    mod_pos = -1;
    mod_num = 1;
    mod_text = "";
    headgroup_decorators->clear();
    add_omega_linoleoyloxy_Cer = false;
}
    
void LipidMapsParserEventHandler::set_molecular_subspecies_level(TreeNode* node){
    set_lipid_level(MOLECULAR_SPECIES);
}

void LipidMapsParserEventHandler::pure_fa(TreeNode* node){
    head_group = "FA";
}

    
    
void LipidMapsParserEventHandler::mediator_event(TreeNode* node){
    use_head_group = true;
    head_group = node->get_text();
}


void LipidMapsParserEventHandler::set_isomeric_level(TreeNode* node){
    db_position = 0;
    db_cistrans = "";
}


const map<string, int> LipidMapsParserEventHandler::acer_heads{
    {"1-O-myristoyl", 14},
    {"1-O-palmitoyl", 16},
    {"1-O-stearoyl", 18},
    {"1-O-eicosanoyl", 20},
    {"1-O-behenoyl", 22},
    {"1-O-lignoceroyl", 24},
    {"1-O-cerotoyl", 26},
    {"1-O-pentacosanoyl", 25},
    {"1-O-carboceroyl", 28},
    {"1-O-tricosanoyl", 30},
    {"1-O-lignoceroyl-omega-linoleoyloxy", 24},
    {"1-O-stearoyl-omega-linoleoyloxy", 18}};
    
    
void LipidMapsParserEventHandler::add_ACer(TreeNode *node){
    string head = node->get_text();
    head_group = "ACer";
    
    if (uncontains_val(acer_heads, head)){
        throw LipidException("ACer head group '" + head + "' unknown");
    }
    
    HeadgroupDecorator *hgd = new HeadgroupDecorator("decorator_acyl", -1, 1, 0, true);
    int acer_num = acer_heads.at(head);
    hgd->functional_groups->insert({"decorator_acyl", vector<FunctionalGroup*>{new FattyAcid("FA", acer_num)}});
    headgroup_decorators->push_back(hgd);
    
    if (head == "1-O-lignoceroyl-omega-linoleoyloxy" || head == "1-O-stearoyl-omega-linoleoyloxy"){
        add_omega_linoleoyloxy_Cer = true;
    }
}
        
        
void LipidMapsParserEventHandler::add_db_position(TreeNode* node){
    if (current_fa != NULL){
        current_fa->double_bonds->double_bond_positions.insert({db_position, db_cistrans});
        if (db_cistrans != "E" && db_cistrans != "Z") set_lipid_level(STRUCTURE_DEFINED);
    }
}


void LipidMapsParserEventHandler::add_db_position_number(TreeNode* node){
    db_position = atoi(node->get_text().c_str());
}


void LipidMapsParserEventHandler::add_cistrans(TreeNode* node){
    db_cistrans = node->get_text();
}
    
    
void LipidMapsParserEventHandler::set_head_group_name(TreeNode* node){
    if (head_group.length() == 0) head_group = node->get_text();
}


void LipidMapsParserEventHandler::set_omega_head_group_name(TreeNode* node){
    add_omega_linoleoyloxy_Cer = true;
    set_head_group_name(node);
}

    
    
void LipidMapsParserEventHandler::set_species_level(TreeNode* node){
    set_lipid_level(SPECIES);
}
   
   
   
void LipidMapsParserEventHandler::set_structural_subspecies_level(TreeNode* node){
    level = min(level, STRUCTURE_DEFINED);
}


void LipidMapsParserEventHandler::set_mod(TreeNode* node){
    mod_text = "";
    mod_pos = -1;
    mod_num = 1;
}


void LipidMapsParserEventHandler::set_mod_text(TreeNode* node){
    mod_text = node->get_text();
}


void LipidMapsParserEventHandler::set_mod_pos(TreeNode* node){
    mod_pos = atoi(node->get_text().c_str());
}


void LipidMapsParserEventHandler::set_mod_num(TreeNode* node){
    mod_num = atoi(node->get_text().c_str());
}   
    
    
void LipidMapsParserEventHandler::add_functional_group(TreeNode* node){
    if (mod_text != "Cp"){
        FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group(mod_text);
        functional_group->position = mod_pos;
        functional_group->count = mod_num;
        string fg_name = functional_group->name;
        if (uncontains_val_p(current_fa->functional_groups, fg_name)) current_fa->functional_groups->insert({fg_name, vector<FunctionalGroup*>()});
        current_fa->functional_groups->at(fg_name).push_back(functional_group);
    }
    else {
        current_fa->num_carbon += 1;
        Cycle *cycle = new Cycle(3, mod_pos, mod_pos + 2);
        if (uncontains_val_p(current_fa->functional_groups, "cy")) current_fa->functional_groups->insert({"cy", vector<FunctionalGroup*>()});
        current_fa->functional_groups->at("cy").push_back(cycle);
    }
}


void LipidMapsParserEventHandler::add_glyco(TreeNode* node){
    string glyco_name = node->get_text();
    HeadgroupDecorator *functional_group = 0;
    try {
        functional_group = (HeadgroupDecorator*)KnownFunctionalGroups::get_functional_group(glyco_name);
    }
    catch (...){
        throw LipidParsingException("Carbohydrate '" + glyco_name + "' unknown");
    }
    
    functional_group->elements->at(ELEMENT_O) -= 1;
    headgroup_decorators->push_back(functional_group);
}

        
        
void LipidMapsParserEventHandler::new_fa(TreeNode *node) {
    db_numbers = -1;
    current_fa = new FattyAcid("FA");
}
    
    

void LipidMapsParserEventHandler::new_lcb(TreeNode *node) {
    lcb = new FattyAcid("LCB");
    lcb->set_type(LCB_REGULAR);
    current_fa = lcb;
    set_lipid_level(STRUCTURE_DEFINED);
}
        
        

void LipidMapsParserEventHandler::clean_lcb(TreeNode *node) {
    if (db_numbers > -1 && db_numbers != current_fa->double_bonds->get_num()){
        throw LipidException("Double bond count does not match with number of double bond positions");
    }
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
    current_fa = NULL;
}
    
    
        

void LipidMapsParserEventHandler::append_fa(TreeNode *node) {
    if (db_numbers > -1 && db_numbers != current_fa->double_bonds->get_num()){
        throw LipidException("Double bond count does not match with number of double bond positions");
    }
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
    
    if (current_fa->num_carbon == 0){
        omit_fa = true;
    }
    fa_list->push_back(current_fa);
    current_fa = NULL;
}
    
    
void LipidMapsParserEventHandler::add_ether(TreeNode* node){
    string ether = node->get_text();
    if (ether == "O-") current_fa->lipid_FA_bond_type = ETHER_PLASMANYL;
    else if (ether == "P-") current_fa->lipid_FA_bond_type = ETHER_PLASMENYL;
}
    
    
    
void LipidMapsParserEventHandler::add_hydroxyl(TreeNode* node){
    int num_h = atoi(node->get_text().c_str());
    
    if (sp_regular_lcb()) num_h -= 1;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("OH").push_back(functional_group);
}


    
void LipidMapsParserEventHandler::add_hydroxyl_lcb(TreeNode* node){
    string hydroxyl = node->get_text();
    int num_h = 0;
    if (hydroxyl == "m") num_h = 1;
    else if (hydroxyl == "d") num_h = 2;
    else if (hydroxyl == "t") num_h = 3;
    
    if (sp_regular_lcb()) num_h -= 1;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("OH").push_back(functional_group);
}
    
    
void LipidMapsParserEventHandler::add_double_bonds(TreeNode* node){
    current_fa->double_bonds->num_double_bonds += atoi(node->get_text().c_str());
}
    
    
void LipidMapsParserEventHandler::add_carbon(TreeNode* node){
    current_fa->num_carbon = atoi(node->get_text().c_str());
}
    
    

void LipidMapsParserEventHandler::build_lipid(TreeNode* node){
    if (omit_fa && head_group_exceptions.find(head_group) != head_group_exceptions.end()){
        head_group = "L" + head_group;
    }
    
    if (lcb != NULL){
        set_lipid_level(STRUCTURE_DEFINED);
        fa_list->insert(fa_list->begin(), lcb);
    }
    
    if (add_omega_linoleoyloxy_Cer){
        if (fa_list->size() != 2){
            throw LipidException("omega-linoleoyloxy-Cer with a different combination to one long chain base and one fatty acyl chain unknown");
        }
        if (uncontains_val_p(fa_list->back()->functional_groups, "acyl")) fa_list->back()->functional_groups->insert({"acyl", vector<FunctionalGroup*>()});
        
        DoubleBonds* db = new DoubleBonds(2);
        db->double_bond_positions.insert({9, "Z"});
        db->double_bond_positions.insert({12, "Z"});
        fa_list->back()->functional_groups->at("acyl").push_back(new AcylAlkylGroup(new FattyAcid("FA", 18, db)));
    }
    
    Headgroup* headgroup = prepare_headgroup_and_checks();
    
    LipidAdduct *lipid = new LipidAdduct();
    lipid->lipid = assemble_lipid(headgroup);
    lipid->adduct = adduct;
    BaseParserEventHandler<LipidAdduct*>::content = lipid;
}
    
    

void LipidMapsParserEventHandler::new_adduct(TreeNode *node) {
    adduct = new Adduct("", "");
}
    
    

void LipidMapsParserEventHandler::add_adduct(TreeNode *node) {
    adduct->adduct_string = node->get_text();
}
    
    

void LipidMapsParserEventHandler::add_charge(TreeNode *node) {
    adduct->charge = atoi(node->get_text().c_str());
}
    
    

void LipidMapsParserEventHandler::add_charge_sign(TreeNode *node) {
    string sign = node->get_text();
    if (sign == "+") adduct->set_charge_sign(1);
    else if (sign == "-") adduct->set_charge_sign(-1);
}
    
        
