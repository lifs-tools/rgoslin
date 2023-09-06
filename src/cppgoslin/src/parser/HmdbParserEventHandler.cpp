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


#include "cppgoslin/parser/HmdbParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<LipidAdduct*>::registered_events->insert({x, bind(&HmdbParserEventHandler::y, this, placeholders::_1)})
    

HmdbParserEventHandler::HmdbParserEventHandler() : LipidBaseParserEventHandler() {
    
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
    reg("st_sub2_hg_pre_event", set_head_group_name);
    reg("ganglioside_names_pre_event", set_head_group_name);
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
    reg("fa_lcb_suffix_type_pre_event", add_one_hydroxyl);
    reg("interlink_fa_pre_event", interlink_fa);
    reg("lipid_suffix_pre_event", lipid_suffix);
    reg("methyl_pre_event", add_methyl);
    reg("furan_fa_pre_event", furan_fa);
    reg("furan_fa_post_event", furan_fa_post);
    reg("furan_fa_mono_pre_event", furan_fa_mono);
    reg("furan_fa_di_pre_event", furan_fa_di);
    reg("furan_first_number_pre_event", furan_fa_first_number);
    reg("furan_second_number_pre_event", furan_fa_second_number);
    
    reg("adduct_info_pre_event", new_adduct);
    reg("adduct_pre_event", add_adduct);
    reg("charge_pre_event", add_charge);
    reg("charge_sign_pre_event", add_charge_sign);
        
    reg("fa_lcb_suffix_types_pre_event", register_suffix_type);
    reg("fa_lcb_suffix_position_pre_event", register_suffix_pos);
    reg("fa_synonym_pre_event", register_fa_synonym);
    
}


HmdbParserEventHandler::~HmdbParserEventHandler(){
}


void HmdbParserEventHandler::reset_lipid(TreeNode *node) {
    level = FULL_STRUCTURE;
    head_group = "";
    lcb = NULL;
    fa_list->clear();
    current_fa = NULL;
    use_head_group = false;
    db_position = 0;
    db_cistrans = "";
    adduct = 0;
    furan.remove_all();
    headgroup_decorators->clear();
    func_type = "";
}


void HmdbParserEventHandler::HmdbParserEventHandler::register_suffix_type(TreeNode* node){
        func_type = node->get_text();
        if (func_type != "me" && func_type != "OH" && func_type != "O"){
            throw LipidException("Unknown functional abbreviation: " + func_type);
        }
        if (func_type == "me") func_type = "Me";
        else if (func_type == "O") func_type = "oxo";
}
    
    
    
void HmdbParserEventHandler::register_suffix_pos(TreeNode* node){
        int pos = node->get_int();
        FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group(func_type);
        functional_group->position = pos;
        if (uncontains_val_p(current_fa->functional_groups, func_type)) current_fa->functional_groups->insert({func_type, vector<FunctionalGroup*>()});
        current_fa->functional_groups->at(func_type).push_back(functional_group);
}



void HmdbParserEventHandler::register_fa_synonym(TreeNode* node){
    string mediator_name = node->get_text();
     
    current_fa = resolve_fa_synonym(mediator_name);
    /*
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
    */
}



void HmdbParserEventHandler::HmdbParserEventHandler::set_isomeric_level(TreeNode* node){
    db_position = 0;
    db_cistrans = "";
}



void HmdbParserEventHandler::add_db_position(TreeNode* node){
    if (current_fa != NULL){
        current_fa->double_bonds->double_bond_positions.insert({db_position, db_cistrans});
        if (db_cistrans != "E" && db_cistrans != "Z") set_lipid_level(STRUCTURE_DEFINED);
    }
}



void HmdbParserEventHandler::add_db_position_number(TreeNode* node){
    db_position = atoi(node->get_text().c_str());
}



void HmdbParserEventHandler::add_cistrans(TreeNode* node){
    db_cistrans = node->get_text();
}



void HmdbParserEventHandler::set_head_group_name(TreeNode *node) {
    head_group = node->get_text();
}


void HmdbParserEventHandler::set_species_level(TreeNode *node) {
    set_lipid_level(SPECIES);
}
    



void HmdbParserEventHandler::set_molecular_level(TreeNode *node) {
    set_lipid_level(MOLECULAR_SPECIES);
}


void HmdbParserEventHandler::mediator_event(TreeNode* node){
    use_head_group = true;
    head_group = node->get_text();
}
    
    

void HmdbParserEventHandler::new_fa(TreeNode *node) {
    current_fa = new FattyAcid("FA");
}
    
    

void HmdbParserEventHandler::new_lcb(TreeNode *node) {
    lcb = new FattyAcid("LCB");
    lcb->set_type(LCB_REGULAR);
    current_fa = lcb;
    set_lipid_level(STRUCTURE_DEFINED);
}
        
        

void HmdbParserEventHandler::clean_lcb(TreeNode *node) {
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }
    current_fa = NULL;
}
    
    
        

void HmdbParserEventHandler::append_fa(TreeNode *node) {
    if (current_fa->double_bonds->get_num() < 0){
        throw LipidException("Double bond count does not match with number of double bond positions");
    }
    if (current_fa->double_bonds->double_bond_positions.size() == 0 && current_fa->double_bonds->get_num() > 0){
        set_lipid_level(SN_POSITION);
    }

    fa_list->push_back(current_fa);
    current_fa = NULL;
}
    
    

void HmdbParserEventHandler::build_lipid(TreeNode *node) {
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
    
    

void HmdbParserEventHandler::add_ether(TreeNode *node) {
    string ether = node->get_text();
    if (ether == "O-" || ether == "o-") current_fa->lipid_FA_bond_type = ETHER_PLASMANYL;
    else if (ether == "P-") current_fa->lipid_FA_bond_type = ETHER_PLASMENYL;
    else throw UnsupportedLipidException("Fatty acyl chain of type '" + ether + "' is currently not supported");
}
    
    

void HmdbParserEventHandler::add_methyl(TreeNode *node) {
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("Me");
    functional_group->position = current_fa->num_carbon - (node->get_text() == "i-" ? 1 : 2);
    current_fa->num_carbon -= 1;
    if (uncontains_val_p(current_fa->functional_groups, "Me")) current_fa->functional_groups->insert({"Me", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("Me").push_back(functional_group);
}
    
    

void HmdbParserEventHandler::add_hydroxyl(TreeNode *node) {
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


void HmdbParserEventHandler::add_one_hydroxyl(TreeNode *node) {
    if (contains_val_p(current_fa->functional_groups, "OH") && current_fa->functional_groups->at("OH").at(0)->position == -1){
        current_fa->functional_groups->at("OH").at(0)->count += 1;
    }
    else {
        FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
        if (uncontains_val_p(current_fa->functional_groups, "OH")) current_fa->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
        current_fa->functional_groups->at("OH").push_back(functional_group);
    }
}
    

void HmdbParserEventHandler::add_double_bonds(TreeNode *node) {
    current_fa->double_bonds->num_double_bonds = atoi(node->get_text().c_str());
}
    
    

void HmdbParserEventHandler::add_carbon(TreeNode *node) {
    current_fa->num_carbon += atoi(node->get_text().c_str());
}
    

void HmdbParserEventHandler::furan_fa(TreeNode *node) {
    furan.remove_all();
}


void HmdbParserEventHandler::furan_fa_post(TreeNode *node) {
    int l = 4 + furan.get_int("len_first") + furan.get_int("len_second");
    current_fa->num_carbon = l;
    
    int start = 1 + furan.get_int("len_first");
    int end = 3 + start;
    DoubleBonds *cyclo_db = new DoubleBonds(2);
    cyclo_db->double_bond_positions.insert({start, "E"});
    cyclo_db->double_bond_positions.insert({2 + start, "E"});
    
    map<string, vector<FunctionalGroup*> > *cyclo_fg = new map<string, vector<FunctionalGroup*> >();
    cyclo_fg->insert({"Me", vector<FunctionalGroup*>()});
    
    if (furan.get_string("type") == "m"){
        FunctionalGroup *fg = KnownFunctionalGroups::get_functional_group("Me");
        fg->position = 1 + start;
        cyclo_fg->at("Me").push_back(fg);
    }
        
    else if (furan.get_string("type") == "d"){
        FunctionalGroup *fg = KnownFunctionalGroups::get_functional_group("Me");
        fg->position = 1 + start;
        cyclo_fg->at("Me").push_back(fg);
        fg = KnownFunctionalGroups::get_functional_group("Me");
        fg->position = 2 + start;
        cyclo_fg->at("Me").push_back(fg);
    }
    
    vector<Element> *bridge_chain = new vector<Element>{ELEMENT_O};
    Cycle *cycle = new Cycle(end - start + 1 + bridge_chain->size(), start, end, cyclo_db, cyclo_fg, bridge_chain);
    current_fa->functional_groups->insert({"cy", vector<FunctionalGroup*>()});
    current_fa->functional_groups->at("cy").push_back(cycle);
}



void HmdbParserEventHandler::furan_fa_mono(TreeNode *node) {
    furan.set_string("type", "m");
}



void HmdbParserEventHandler::furan_fa_di(TreeNode *node) {
    furan.set_string("type", "d");
}



void HmdbParserEventHandler::furan_fa_first_number(TreeNode *node) {
    furan.set_int("len_first", atoi(node->get_text().c_str()));
}



void HmdbParserEventHandler::furan_fa_second_number(TreeNode *node) {
    furan.set_int("len_second", atoi(node->get_text().c_str()));
    
}
    

void HmdbParserEventHandler::interlink_fa(TreeNode *node) {
    throw UnsupportedLipidException("Interconnected fatty acyl chains are currently not supported");
}
    

void HmdbParserEventHandler::lipid_suffix(TreeNode *node) {
    //throw UnsupportedLipidException("Lipids with suffix '" + node->get_text() + "' are currently not supported");
}
    
    

void HmdbParserEventHandler::new_adduct(TreeNode *node) {
    if (!adduct) adduct = new Adduct("", "");
}
    
    

void HmdbParserEventHandler::add_adduct(TreeNode *node) {
    adduct->adduct_string = node->get_text();
}
    
    

void HmdbParserEventHandler::add_charge(TreeNode *node) {
    adduct->charge = atoi(node->get_text().c_str());
}
    
    

void HmdbParserEventHandler::add_charge_sign(TreeNode *node) {
    string sign = node->get_text();
    if (sign == "+") adduct->set_charge_sign(1);
    else if (sign == "-") adduct->set_charge_sign(-1);
    if (adduct->charge == 0) adduct->charge = 1;
}

        

        
