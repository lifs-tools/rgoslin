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
const map<string, int> GoslinParserEventHandler::mediator_trivial{{"Palmitic acid", 0}, {"Linoleic acid", 1}, {"AA", 2}, {"ALA", 3}, {"EPA", 4}, {"DHA", 5}, {"LTB4", 6}, {"Resolvin D3", 7}, {"Maresin 1", 8},  {"Resolvin D2", 9}, {"Resolvin D5", 10}, {"Resolvin D1", 11}, {"TXB1", 12}, {"TXB2", 13}, {"TXB3", 14}, {"PGF2alpha", 15}, {"PGD2", 16}, {"PGE2", 17}, {"PGB2", 18}, {"15d-PGJ2", 19}};

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
    reg("mediator_position_pre_event", set_mediator_function_position);
    reg("mediator_functional_group_post_event", add_mediator_function);
    reg("mediator_suffix_pre_event", add_mediator_suffix);
    reg("mediator_tetranor_pre_event", set_mediator_tetranor);
        
    debug = "";
}


GoslinParserEventHandler::~GoslinParserEventHandler(){
}


void GoslinParserEventHandler::set_mediator(TreeNode *node){
    head_group = "FA";
    current_fa = new FattyAcid("FA");
    fa_list->push_back(current_fa);
    set_lipid_level(STRUCTURE_DEFINED);
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
    current_fa->num_carbon += GoslinParserEventHandler::mediator_FA.at(node->get_text());
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
    FunctionalGroup* functional_group = 0;
    string fg = "";
    if (mediator_function == "H"){
        functional_group = KnownFunctionalGroups::get_functional_group("OH");
        fg = "OH";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "Oxo"){
        functional_group = KnownFunctionalGroups::get_functional_group("oxo");
        fg = "oxo";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "Hp"){
        functional_group = KnownFunctionalGroups::get_functional_group("OOH");
        fg = "OOH";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "E" || mediator_function == "Ep"){
        functional_group = KnownFunctionalGroups::get_functional_group("Ep");
        fg = "Ep";
        if (mediator_function_positions.size() > 0) functional_group->position = mediator_function_positions[0];
    }
        
    else if (mediator_function == "DH" || mediator_function == "DiH" || mediator_function == "diH"){
        functional_group = KnownFunctionalGroups::get_functional_group("OH");
        fg = "OH";
        if (mediator_function_positions.size() > 0){
            functional_group->position = mediator_function_positions[0];
            FunctionalGroup* functional_group2 = KnownFunctionalGroups::get_functional_group("OH");
            functional_group2->position = mediator_function_positions[1];
            current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
            current_fa->functional_groups->at("OH").push_back(functional_group2);
        }
    }
        
    if (uncontains_val_p(current_fa->functional_groups, fg)) current_fa->functional_groups->insert({fg, vector<FunctionalGroup*>()});
    current_fa->functional_groups->at(fg).push_back(functional_group);
}

void GoslinParserEventHandler::set_trivial_mediator(TreeNode *node){
    head_group = "FA";
    string mediator_name = node->get_text();
     
    switch(GoslinParserEventHandler::mediator_trivial.at(mediator_name)){
        case 0: // Palmitic acid
            current_fa = new FattyAcid("FA", 16);
            break;
            
        case 1: // Linoleic acid":
            current_fa = new FattyAcid("FA", 18, new DoubleBonds(2));
            break;
            
        case 2: // AA":
            current_fa = new FattyAcid("FA", 20, new DoubleBonds(4));
            break;
            
        case 3: // ALA":
            current_fa = new FattyAcid("FA", 18, new DoubleBonds(3));
            break;
            
        case 4: // EPA":
            current_fa = new FattyAcid("FA", 20, new DoubleBonds(5));
            break;
            
        case 5: // DHA":
            current_fa = new FattyAcid("FA", 22, new DoubleBonds(6));
            break;
            
        case 6: // LTB4
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 5;
                f2->position = 12;
                current_fa = new FattyAcid("FA", 20, new DoubleBonds({{6, "Z"}, {8, "E"}, {10, "E"}, {14, "Z"}}), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
            }
            break;
            
        case 7: // Resolvin D3
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 4;
                f2->position = 11;
                f3->position = 17;
                current_fa = new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2, f3}}});
            }
            break;
            
        case 8: // Maresin 1
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 4;
                f2->position = 14;
                current_fa = new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
            }
            break;
            
        case 9: // Resolvin D2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 4;
                f2->position = 16;
                f3->position = 17;
                current_fa = new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2, f3}}});
            }
            break;
            
        case 10: // Resolvin D5
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 7;
                f2->position = 17;
                current_fa = new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
            }
            break;
            
        case 11: // Resolvin D1
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 7;
                f2->position = 8;
                f3->position = 17;
                current_fa = new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2, f3}}});
            }
            break;
            
        case 12: // TXB1
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f4 = KnownFunctionalGroups::get_functional_group("oxy");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                f4->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}, {"oxy", {f4}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 13: // TXB2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f4 = KnownFunctionalGroups::get_functional_group("oxy");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                f4->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}, {"oxy", {f4}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 14: // TXB3
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f4 = KnownFunctionalGroups::get_functional_group("oxy");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                f4->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}, {"oxy", {f4}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(3), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 15: // PGF2alpha
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 16: // PGD2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("oxo");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2}}, {"oxo", {f3}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 17: // PGE2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f3}}, {"oxy", {f2}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 18: // PGB2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 15;
                f2->position = 9;
                Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"OH", {f2}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 19: // 15d-PGJ2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
                f1->position = 15;
                f2->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"oxo", {f2}}});
                current_fa = new FattyAcid("FA", 20, new DoubleBonds(3), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
    }
    
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
    current_fa->double_bonds->num_double_bonds = atoi(node->get_text().c_str());
}
    
    

void GoslinParserEventHandler::add_carbon(TreeNode *node) {
    current_fa->num_carbon = atoi(node->get_text().c_str());
}
    
    

void GoslinParserEventHandler::add_hydroxyl(TreeNode *node) {
    int num_h = atoi(node->get_text().c_str());
    if (sp_regular_lcb()) num_h -= 1;
    if (num_h <= 0) return;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("OH").push_back(functional_group);
    level = min(level, STRUCTURE_DEFINED);
}
    
    

void GoslinParserEventHandler::new_adduct(TreeNode *node) {
    adduct = new Adduct("", "");
}
    
    

void GoslinParserEventHandler::add_adduct(TreeNode *node) {
    adduct->adduct_string = node->get_text();
}
    
    

void GoslinParserEventHandler::add_charge(TreeNode *node) {
    adduct->charge = atoi(node->get_text().c_str());
}
    
    

void GoslinParserEventHandler::add_charge_sign(TreeNode *node) {
    string sign = node->get_text();
    if (sign == "+") adduct->set_charge_sign(1);
    else if (sign == "-") adduct->set_charge_sign(-1);
}
        

        
