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


#include "cppgoslin/parser/GoslinParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<LipidAdduct*>::registered_events->insert({x, bind(&GoslinParserEventHandler::y, this, placeholders::_1)})
   
   
const map<string, int> GoslinParserEventHandler::mediator_FA{{"H", 17}, {"O", 18}, {"E", 20}, {"Do", 22}};
const map<string, int> GoslinParserEventHandler::mediator_DB{{"M", 1}, {"D", 2}, {"Tr", 3}, {"T", 4}, {"P", 5}, {"H", 6}};


GoslinParserEventHandler::GoslinParserEventHandler() : LipidBaseParserEventHandler() {    
    reg("lipid_pre_event", reset_lipid);
    reg("lipid_post_event", build_lipid);
    reg("hg_cl_pre_event", set_head_group_name);
    reg("hg_mlcl_pre_event", set_head_group_name);
    reg("hg_pl_pre_event", set_head_group_name);
    reg("hg_lpl_pre_event", set_head_group_name);
    reg("hg_lsl_pre_event", set_head_group_name);
    reg("hg_so_lsl_pre_event", set_head_group_name);
    reg("hg_dsl_pre_event", set_head_group_name);
    reg("st_pre_event", set_head_group_name);
    reg("hg_ste_pre_event", set_head_group_name);
    reg("hg_stes_pre_event", set_head_group_name);
    reg("hg_mgl_pre_event", set_head_group_name);
    reg("hg_dgl_pre_event", set_head_group_name);
    reg("hg_sgl_pre_event", set_head_group_name);
    reg("hg_tgl_pre_event", set_head_group_name);
    reg("hg_dlcl_pre_event", set_head_group_name);
    reg("hg_sac_di_pre_event", set_head_group_name);
    reg("hg_sac_f_pre_event", set_head_group_name);
    reg("hg_tpl_pre_event", set_head_group_name);  
    reg("gl_species_pre_event", set_species_level);
    reg("pl_species_pre_event", set_species_level);
    reg("sl_species_pre_event", set_species_level);
    reg("fa2_unsorted_pre_event", set_molecular_subspecies_level);
    reg("fa3_unsorted_pre_event", set_molecular_subspecies_level);
    reg("fa4_unsorted_pre_event", set_molecular_subspecies_level);
    reg("slbpa_pre_event", set_molecular_subspecies_level);
    reg("dlcl_pre_event", set_molecular_subspecies_level);
    reg("mlcl_pre_event", set_molecular_subspecies_level);
    reg("lcb_pre_event", new_lcb);
    reg("lcb_post_event", clean_lcb);
    reg("fa_pre_event", new_fa);
    reg("fa_post_event", append_fa);
    reg("db_single_position_pre_event", set_isomeric_level);
    reg("db_single_position_post_event", add_db_position);
    reg("db_position_number_pre_event", add_db_position_number);
    reg("cistrans_pre_event", add_cistrans);
    reg("ether_pre_event", add_ether);
    reg("old_hydroxyl_pre_event", add_old_hydroxyl);
    reg("db_count_pre_event", add_double_bonds);
    reg("carbon_pre_event", add_carbon);
    reg("hydroxyl_pre_event", add_hydroxyl);
    reg("adduct_info_pre_event", new_adduct);
    reg("adduct_pre_event", add_adduct);
    reg("charge_pre_event", add_charge);
    reg("charge_sign_pre_event", add_charge_sign);
    reg("lpl_pre_event", set_molecular_subspecies_level);
    reg("plasmalogen_pre_event", set_plasmalogen);
    reg("mediator_pre_event", set_mediator);
    reg("mediator_post_event", add_mediator);
    reg("unstructured_mediator_pre_event", set_unstructured_mediator);
    reg("trivial_mediator_pre_event", set_trivial_mediator);
    reg("mediator_carbon_pre_event", set_mediator_carbon);
    reg("mediator_db_pre_event", set_mediator_db);
    reg("mediator_mono_functions_pre_event", set_mediator_function);
    reg("mediator_di_functions_pre_event", set_mediator_function);
    reg("mediator_tri_functions_pre_event", set_mediator_function);
    reg("mediator_position_pre_event", set_mediator_function_position);
    reg("mediator_functional_group_post_event", add_mediator_function);
    reg("mediator_suffix_pre_event", add_mediator_suffix);
    reg("mediator_tetranor_pre_event", set_mediator_tetranor);
    reg("isotope_pair_pre_event", new_adduct);
    reg("isotope_element_pre_event", set_heavy_d_element);
    reg("isotope_number_pre_event", set_heavy_d_number);
    reg("heavy_pre_event", new_adduct);
    reg("adduct_heavy_element_pre_event", set_heavy_element);
    reg("adduct_heavy_number_pre_event", set_heavy_number);
    reg("adduct_heavy_component_post_event", add_heavy_component);
    
    reg("prostaglandin_number_pre_event", set_prostaglandin_number);
    reg("prostaglandin_type_pre_event", set_prostaglandin_type);
    reg("prostaglandin_post_event", add_prostaglandin);
        
    debug = "";
}


GoslinParserEventHandler::~GoslinParserEventHandler(){
}


void GoslinParserEventHandler::reset_lipid(TreeNode *node) {
    level = FULL_STRUCTURE;
    head_group = "";
    lcb = NULL;
    fa_list->clear();
    current_fa = NULL;
    adduct = NULL;
    db_position = 0;
    db_cistrans = "";
    unspecified_ether = false;
    plasmalogen = 0;
    mediator_function = "";
    mediator_function_positions.clear();
    mediator_suffix = false;
    use_head_group = false;
    headgroup_decorators->clear();
    heavy_element = ELEMENT_C;
    heavy_element_number = 0;
    trivial_mediator = false;
    prostaglandin_type = "";
    prostaglandin_number = "";
}



        
        
        
void GoslinParserEventHandler::add_prostaglandin(TreeNode *node){
    set<string> pg_types = {"B", "D", "E", "F", "J", "K"};
    set<string> pg_numbers = {"1", "2", "3"};
    if (uncontains_val(pg_types, prostaglandin_type) || uncontains_val(pg_numbers, prostaglandin_number)) return;
    
    DoubleBonds *db = nullptr;
    FattyAcid *tmp_fa = current_fa;
    
    if (prostaglandin_number == "1") db = new DoubleBonds({{13, "E"}});
    else if (prostaglandin_number == "2") db = new DoubleBonds({{5, "Z"}, {13, "E"}});
    else if (prostaglandin_number == "3") db = new DoubleBonds({{5, "Z"}, {13, "E"}, {17, "Z"}});
        
    if (prostaglandin_type == "B"){
        FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
        f1->position = 15;
        f2->position = 9;
        Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds({{8, ""}}), new map<string, vector<FunctionalGroup*>>{{"OH", {f2}}});
        current_fa = new FattyAcid("FA", 20, db, new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
    }
    else if (prostaglandin_type == "D"){
        FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("oxo");
        f1->position = 15;
        f2->position = 9;
        f3->position = 11;
        Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2}}, {"oxo", {f3}}});
        current_fa = new FattyAcid("FA", 20, db, new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
    }
    else if (prostaglandin_type == "E"){
        FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
        FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
        f1->position = 15;
        f2->position = 9;
        f3->position = 11;
        Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f3}}, {"oxy", {f2}}});
        current_fa = new FattyAcid("FA", 20, db, new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
    }
    else if (prostaglandin_type == "F"){
        FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
        f1->position = 15;
        f2->position = 9;
        f3->position = 11;
        Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}});
        current_fa = new FattyAcid("FA", 20, db, new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
    }
    else if (prostaglandin_type == "J"){
        FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
        f1->position = 15;
        f2->position = 11;
        Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds({{9, ""}}), new map<string, vector<FunctionalGroup*>>{{"oxo", {f2}}});
        current_fa = new FattyAcid("FA", 20, db, new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
    }
    else if (prostaglandin_type == "K"){
        FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
        FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
        FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("oxo");
        f1->position = 15;
        f2->position = 9;
        f3->position = 11;
        Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"oxo", {f2, f3}}});
        current_fa = new FattyAcid("FA", 20, db, new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
    }
    else {
        delete db;
        return;
    }
    
    clean_mediator(tmp_fa);
    fa_list->clear();
    fa_list->push_back(current_fa);
    mediator_suffix = true;
    prostaglandin_type = "";
    prostaglandin_number = "";
    delete tmp_fa;
}


bool recursive_deletion(string fg_name, FunctionalGroup *fg, map<string, vector<FunctionalGroup*>> *functional_groups){
    vector<string> del_functions;
    bool return_value = false;
    for (auto &kv : *functional_groups){
        vector<int> del_fg;
        for (int i = 0; i < (int)kv.second.size(); ++i){
            if (kv.second.at(i)->position == fg->position) del_fg.push_back(i);
        }
        if (!del_fg.empty()){
            for (int i = (int)del_fg.size() - 1; i >= 0; --i){
                delete kv.second.at(i);
                kv.second.erase(kv.second.begin() + i);
                if (kv.second.empty()) del_functions.push_back(kv.first);
            }
            if (uncontains_val_p(functional_groups, fg_name)) functional_groups->insert({fg_name, vector<FunctionalGroup*>()});
            functional_groups->at(fg_name).push_back(fg);
            return_value = true;
        }
        for (auto fg_curr : kv.second){
            if (return_value) break;
            if (recursive_deletion(fg_name, fg, fg_curr->functional_groups)) return_value = true;
        }
    }
    for (auto del_fg_name : del_functions) functional_groups->erase(del_fg_name);
    return return_value;
}
        
        
void remove_deoxy(map<string, vector<FunctionalGroup*> > *functional_groups){
    if (contains_val_p(functional_groups, "d")){
        for (auto fg : functional_groups->at("d")) delete fg;
        functional_groups->erase("d");
    }
    for (auto &kv : *functional_groups){
        for (auto fg : kv.second)
            remove_deoxy(fg->functional_groups);
    }
}


void GoslinParserEventHandler::clean_mediator(FattyAcid *fa){
    if (!fa->functional_groups->empty()){
        for (auto &kv : *(fa->functional_groups)){
            string fg_name = kv.first;
            for (auto fg : kv.second){
                
                if (fg->position >= -1){
                    // erase prevously added functional group, if findable
                    if (!recursive_deletion(fg_name, fg, current_fa->functional_groups)){
                        if (uncontains_val_p(current_fa->functional_groups,  fg_name)) current_fa->functional_groups->insert({fg_name, vector<FunctionalGroup*>()});
                        current_fa->functional_groups->at(fg_name).push_back(fg);
                    }
                }
                else {
                    if (uncontains_val_p(current_fa->functional_groups,  fg_name)) current_fa->functional_groups->insert({fg_name, vector<FunctionalGroup*>()});
                    current_fa->functional_groups->at(fg_name).push_back(fg);
                }
            }
        }
        fa->functional_groups->clear();
    }
    remove_deoxy(current_fa->functional_groups);
    
}

        
void GoslinParserEventHandler::set_prostaglandin_type(TreeNode *node){
    prostaglandin_type = node->get_text();
}
    
    
    
void GoslinParserEventHandler::set_prostaglandin_number(TreeNode *node){
    prostaglandin_number = node->get_text();
}
        
        
        
void GoslinParserEventHandler::set_mediator(TreeNode *node){
    head_group = "FA";
    current_fa = new FattyAcid("FA");
    fa_list->push_back(current_fa);
}
    
    
void GoslinParserEventHandler::set_unstructured_mediator(TreeNode *node){
    head_group = node->get_text();
    use_head_group = true;
    fa_list->clear();
}
        
        
void GoslinParserEventHandler::set_mediator_tetranor(TreeNode *node){
    current_fa->num_carbon -= 4;
}
    

void GoslinParserEventHandler::set_mediator_carbon(TreeNode *node){
    trivial_mediator = true;
    current_fa->num_carbon += GoslinParserEventHandler::mediator_FA.at(node->get_text());
}
        
        
        
void GoslinParserEventHandler::set_heavy_d_element(TreeNode *node){
    adduct->heavy_elements[ELEMENT_H2] = 1;
}
        
        
        
void GoslinParserEventHandler::set_heavy_d_number(TreeNode *node){
    adduct->heavy_elements[ELEMENT_H2] = node->get_int();
}

    

void GoslinParserEventHandler::set_mediator_db(TreeNode *node){
    current_fa->double_bonds->num_double_bonds = GoslinParserEventHandler::mediator_DB.at(node->get_text());
}
    
    
void GoslinParserEventHandler::set_mediator_function(TreeNode *node){
    mediator_function = node->get_text();
}
    
    
void GoslinParserEventHandler::set_mediator_function_position(TreeNode *node){
    mediator_function_positions.push_back(node->get_int());
}
    
    
void GoslinParserEventHandler::add_mediator_function(TreeNode *node){
    set<string> oxos = {"oxo", "Oxo", "OXO", "k", "keto", "K"};
    set<string> dihydro = {"dh", "DH", "DiH", "diH", "dihydro"};
    set<string> trihydro = {"triH", "TriH", "trihydro"};
    
    FunctionalGroup* functional_group = 0;
    string fg = "";
    if (mediator_function == "H"){
        functional_group = KnownFunctionalGroups::get_functional_group("OH");
        fg = "OH";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
    
    else if (mediator_function == "d" || mediator_function == "deoxy"){
        if (mediator_function_positions.size() > 0){
            functional_group = KnownFunctionalGroups::get_functional_group("d");
            functional_group->position = mediator_function_positions[0];
        }
        fg = "d";
    }
        
    else if (contains_val(oxos, mediator_function)){
        functional_group = KnownFunctionalGroups::get_functional_group("oxo");
        fg = "oxo";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "Hp" || mediator_function == "HP"){
        functional_group = KnownFunctionalGroups::get_functional_group("OOH");
        fg = "OOH";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "E" || mediator_function == "Ep"){
        functional_group = KnownFunctionalGroups::get_functional_group("Ep");
        fg = "Ep";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "NO2"){
        functional_group = KnownFunctionalGroups::get_functional_group("NO2");
        fg = "NO2";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (contains_val(dihydro, mediator_function)){
        functional_group = KnownFunctionalGroups::get_functional_group("OH");
        fg = "OH";
        if (mediator_function_positions.size() > 1){
            functional_group->position = mediator_function_positions[0];
            FunctionalGroup* functional_group2 = KnownFunctionalGroups::get_functional_group("OH");
            functional_group2->position = mediator_function_positions[1];
            if (uncontains_val_p(current_fa->functional_groups, fg)) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
            current_fa->functional_groups->at("OH").push_back(functional_group2);
        }
        else {
            functional_group->count = 2;
        }
    }
        
    else if (contains_val(trihydro, mediator_function)){
        functional_group = KnownFunctionalGroups::get_functional_group("OH");
        fg = "OH";
        if (mediator_function_positions.size() > 2){
            functional_group->position = mediator_function_positions[0];
            FunctionalGroup* functional_group2 = KnownFunctionalGroups::get_functional_group("OH");
            functional_group2->position = mediator_function_positions[1];
            FunctionalGroup* functional_group3 = KnownFunctionalGroups::get_functional_group("OH");
            functional_group2->position = mediator_function_positions[2];
            if (uncontains_val_p(current_fa->functional_groups, fg)) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
            current_fa->functional_groups->at("OH").push_back(functional_group2);
            current_fa->functional_groups->at("OH").push_back(functional_group3);
        }
        else {
            functional_group->count = 3;
        }
    }
    
    if (functional_group == nullptr) return;
        
    if (uncontains_val_p(current_fa->functional_groups, fg)) current_fa->functional_groups->insert({fg, vector<FunctionalGroup*>()});
    current_fa->functional_groups->at(fg).push_back(functional_group);
    mediator_function_positions.clear();
    mediator_function = "";
}



void GoslinParserEventHandler::set_trivial_mediator(TreeNode *node){
    head_group = "FA";
    string mediator_name = node->get_text();
    
    FattyAcid *tmp_fa = current_fa;
    current_fa = resolve_fa_synonym(mediator_name);
    clean_mediator(tmp_fa);
    delete tmp_fa;
    
    fa_list->clear();
    fa_list->push_back(current_fa);
    mediator_suffix = true;
}
    
    
void GoslinParserEventHandler::add_mediator_suffix(TreeNode *node){
    mediator_suffix = true;
}
    
    
void GoslinParserEventHandler::add_mediator(TreeNode *node){
    if (!mediator_suffix){
        current_fa->double_bonds->num_double_bonds -= 1;
    }
}


void GoslinParserEventHandler::set_plasmalogen(TreeNode *node) {
    plasmalogen = toupper(node->get_text()[0]);
}


void GoslinParserEventHandler::set_unspecified_ether(TreeNode *node) {
    //unspecified_ether = true;
}

void GoslinParserEventHandler::set_head_group_name(TreeNode *node) {
    head_group = node->get_text();
}


void GoslinParserEventHandler::set_species_level(TreeNode *node) {
    set_lipid_level(SPECIES);
}


void GoslinParserEventHandler::set_isomeric_level(TreeNode* node){
    db_position = 0;
    db_cistrans = "";
}


void GoslinParserEventHandler::add_db_position(TreeNode* node){
    if (current_fa != NULL){
        current_fa->double_bonds->double_bond_positions.insert({db_position, db_cistrans});
        if (db_cistrans != "E" && db_cistrans != "Z") set_lipid_level(STRUCTURE_DEFINED);
    }
}


void GoslinParserEventHandler::add_db_position_number(TreeNode* node){
    db_position = atoi(node->get_text().c_str());
}


void GoslinParserEventHandler::add_cistrans(TreeNode* node){
    db_cistrans = node->get_text();
}
    
    

void GoslinParserEventHandler::set_molecular_subspecies_level(TreeNode *node) {
    set_lipid_level(MOLECULAR_SPECIES);
}
    
    
void GoslinParserEventHandler::new_fa(TreeNode *node) {
    LipidFaBondType lipid_FA_bond_type = ESTER;
    if (unspecified_ether){
        unspecified_ether = false;
        lipid_FA_bond_type = ETHER_UNSPECIFIED;
    }
    current_fa = new FattyAcid("FA", 2, 0, 0, lipid_FA_bond_type);
}
    
    

void GoslinParserEventHandler::new_lcb(TreeNode *node) {
    lcb = new FattyAcid("LCB");
    current_fa = lcb;
    set_lipid_level(STRUCTURE_DEFINED);
    lcb->set_type(LCB_REGULAR);
}
        
        

void GoslinParserEventHandler::clean_lcb(TreeNode *node) {
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
    current_fa = NULL;
}
    
    
        

void GoslinParserEventHandler::append_fa(TreeNode *node) {
    if (current_fa->lipid_FA_bond_type == ETHER_UNSPECIFIED){
        throw LipidException("Lipid with unspecified ether bond cannot be treated properly.");
    }
    
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
        
    if (current_fa->double_bonds->get_num() < 0){
        throw LipidException("Double bond count does not match with number of double bond positions");
    }

    remove_deoxy(current_fa->functional_groups);
    fa_list->push_back(current_fa);
    current_fa = NULL;
    
    if (head_group == "Sa" || head_group == "So" || head_group == "S1P" || head_group == "Sa1P"){
        FattyAcid* fa = fa_list->at(0);
        
        FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
        if (head_group == "Sa" || head_group == "So"){
            functional_group->count = 2;
            fa->lipid_FA_bond_type = LCB_EXCEPTION;
        }
        else {
            functional_group->count = 1;
            fa->lipid_FA_bond_type = LCB_REGULAR;
        }
        if (uncontains_val_p(fa->functional_groups, "OH")) fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
        fa->functional_groups->at("OH").push_back(functional_group);
    }
}
    
    

void GoslinParserEventHandler::build_lipid(TreeNode *node) {
    if (lcb){
        set_lipid_level(STRUCTURE_DEFINED);
        fa_list->insert(fa_list->begin(), lcb);
    }
    
    if (plasmalogen && fa_list->size() && lcb == 0){
        fa_list->at(0)->lipid_FA_bond_type = plasmalogen == 'O' ? ETHER_PLASMANYL : ETHER_PLASMENYL;
    }
    
    Headgroup *headgroup = prepare_headgroup_and_checks();
    
    
    string lipid_name = to_lower(strip(node->get_text(), 1));
    map<string, vector<int>> &trivial_db = TrivialMediators::get_instance().trivial_mediators;
    
    if (trivial_mediator && contains_val(trivial_db, lipid_name)){
        vector<int> &db_pos = trivial_db[lipid_name];
        fa_list->at(0)->double_bonds->num_double_bonds = db_pos.size();
        fa_list->at(0)->double_bonds->double_bond_positions.clear();
        for (auto p : db_pos) fa_list->at(0)->double_bonds->double_bond_positions.insert({p, ""});
        level = FULL_STRUCTURE;
    }
    
    LipidAdduct *lipid = new LipidAdduct();
    lipid->lipid = assemble_lipid(headgroup);
    lipid->adduct = adduct;
    BaseParserEventHandler<LipidAdduct*>::content = lipid;
    
}
    
    

void GoslinParserEventHandler::add_ether(TreeNode *node) {
    string ether = node->get_text();
    if (ether == "a") current_fa->lipid_FA_bond_type = ETHER_PLASMANYL;
    else if (ether == "p"){
        current_fa->lipid_FA_bond_type = ETHER_PLASMENYL;
        current_fa->double_bonds->num_double_bonds = max(0, current_fa->double_bonds->num_double_bonds - 1);
    }
    plasmalogen = 0;
}
    
    

void GoslinParserEventHandler::add_old_hydroxyl(TreeNode *node) {
    string old_hydroxyl = node->get_text();
    int num_h = 0;
    if (old_hydroxyl == "d") num_h = 2;
    else if (old_hydroxyl == "t") num_h = 3;

    if (sp_regular_lcb()) num_h -= 1;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("OH").push_back(functional_group);
}
    
    

void GoslinParserEventHandler::add_double_bonds(TreeNode *node) {
    current_fa->double_bonds->num_double_bonds = node->get_int();
}
    
    

void GoslinParserEventHandler::add_carbon(TreeNode *node) {
    current_fa->num_carbon = node->get_int();
}
    
    

void GoslinParserEventHandler::add_hydroxyl(TreeNode *node) {
    int num_h = node->get_int();
    if (sp_regular_lcb()) num_h -= 1;
    if (num_h <= 0) return;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("OH").push_back(functional_group);
    level = min(level, STRUCTURE_DEFINED);
}
    
    

void GoslinParserEventHandler::new_adduct(TreeNode *node) {
    if (!adduct) adduct = new Adduct("", "");
}
    
    

void GoslinParserEventHandler::add_adduct(TreeNode *node) {
    adduct->adduct_string = node->get_text();
}
    
    

void GoslinParserEventHandler::add_charge(TreeNode *node) {
    adduct->charge = node->get_int();
}
    
    

void GoslinParserEventHandler::add_charge_sign(TreeNode *node) {
    string sign = node->get_text();
    if (sign == "+") adduct->set_charge_sign(1);
    else if (sign == "-") adduct->set_charge_sign(-1);
    if (adduct->charge == 0) adduct->charge = 1;
}
        
        
        
void GoslinParserEventHandler::set_heavy_element(TreeNode *node) {
    heavy_element = heavy_element_table.at(node->get_text());
    heavy_element_number = 1;
}
        
        
void GoslinParserEventHandler::set_heavy_number(TreeNode *node) {
    heavy_element_number = node->get_int();
}
        
        
void GoslinParserEventHandler::add_heavy_component(TreeNode *node) {
    adduct->heavy_elements[heavy_element] += heavy_element_number;
}
        

        
