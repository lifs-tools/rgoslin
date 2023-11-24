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


#include "cppgoslin/parser/ShorthandParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<LipidAdduct*>::registered_events->insert({x, bind(&ShorthandParserEventHandler::y, this, placeholders::_1)})
#define FA_I ("fa" + std::to_string(current_fas.size()))



ShorthandParserEventHandler::ShorthandParserEventHandler() : LipidBaseParserEventHandler() {
    
    reg("lipid_pre_event", reset_lipid);
    reg("lipid_post_event", build_lipid);
    
    // set adduct events
    reg("adduct_info_pre_event", new_adduct);
    reg("adduct_pre_event", add_adduct);
    reg("charge_pre_event", add_charge);
    reg("charge_sign_pre_event", add_charge_sign);
    
    // set species events
    reg("med_species_pre_event", set_species_level);
    reg("gl_species_pre_event", set_species_level);
    reg("gl_molecular_species_pre_event", set_molecular_level);
    reg("pl_species_pre_event", set_species_level);
    reg("pl_molecular_species_pre_event", set_molecular_level);
    reg("sl_species_pre_event", set_species_level);
    reg("pl_single_pre_event", set_molecular_level);
    reg("unsorted_fa_separator_pre_event", set_molecular_level);
    reg("ether_num_pre_event", set_ether_num);
    
    // set head groups events
    reg("med_hg_single_pre_event", set_headgroup_name);
    reg("med_hg_double_pre_event", set_headgroup_name);
    reg("med_hg_triple_pre_event", set_headgroup_name);
    reg("gl_hg_single_pre_event", set_headgroup_name);
    reg("gl_hg_double_pre_event", set_headgroup_name);
    reg("gl_hg_glycosyl_single_pre_event", set_headgroup_name);
    reg("gl_hg_glycosyl_double_pre_event", set_headgroup_name);
    reg("gl_hg_triple_pre_event", set_headgroup_name);
    reg("pl_hg_single_pre_event", set_headgroup_name);
    reg("pl_hg_double_pre_event", set_headgroup_name);
    reg("pl_cpa_head_pre_event", set_headgroup_name);
    reg("pl_hg_quadro_pre_event", set_headgroup_name);
    reg("sl_hg_single_pre_event", set_headgroup_name);
    reg("sl_hg_glyco_pre_event", set_headgroup_name);
    reg("pl_hg_double_fa_hg_pre_event", set_headgroup_name);
    reg("sl_hg_double_name_pre_event", set_headgroup_name);
    reg("st_hg_pre_event", set_headgroup_name);
    reg("st_hg_ester_pre_event", set_headgroup_name);
    reg("hg_pip_pure_m_pre_event", set_headgroup_name);
    reg("hg_pip_pure_d_pre_event", set_headgroup_name);
    reg("hg_pip_pure_t_pre_event", set_headgroup_name);
    reg("hg_PE_PS_pre_event", set_headgroup_name);
    reg("glyco_sphingo_lipid_pre_event", set_glyco_sphingo_lipid);
    reg("carbohydrate_number_pre_event", set_carbohydrate_number);

    // set head group headgroup_decorators
    reg("carbohydrate_sn_pre_event", set_carbohydrate);
    reg("carbohydrate_iso_pre_event", set_carbohydrate);
    reg("carbohydrate_sn_position_pre_event", set_carbohydrate_sn_position);
    reg("carbohydrate_isomeric_pre_event", set_carbohydrate_isomeric);
    
    // fatty acyl events
    reg("lcb_pre_event", new_lcb);
    reg("lcb_post_event", add_fatty_acyl_chain);
    reg("fatty_acyl_chain_pre_event", new_fatty_acyl_chain);
    reg("fatty_acyl_chain_post_event", add_fatty_acyl_chain);
    reg("carbon_pre_event", set_carbon);
    reg("db_count_pre_event", set_double_bond_count);
    reg("db_position_number_pre_event", set_double_bond_position);
    reg("db_single_position_pre_event", set_double_bond_information);
    reg("db_single_position_post_event", add_double_bond_information);
    reg("cistrans_pre_event", set_cistrans);
    reg("ether_type_pre_event", set_ether_type);
    reg("stereo_type_fa_pre_event", set_fatty_acyl_stereo);
    
    // set functional group events
    reg("func_group_data_pre_event", set_functional_group);
    reg("func_group_data_post_event", add_functional_group);
    reg("func_group_pos_number_pre_event", set_functional_group_position);
    reg("func_group_name_pre_event", set_functional_group_name);
    reg("func_group_count_pre_event", set_functional_group_count);
    reg("stereo_type_fg_pre_event", set_functional_group_stereo);
    reg("molecular_func_group_name_pre_event", set_sn_position_func_group);
    reg("fa_db_only_post_event", add_dihydroxyl);
    
    // set cycle events
    reg("func_group_cycle_pre_event", set_cycle);
    reg("func_group_cycle_post_event", add_cycle);
    reg("cycle_start_pre_event", set_cycle_start);
    reg("cycle_end_pre_event", set_cycle_end);
    reg("cycle_number_pre_event", set_cycle_number);
    reg("cycle_db_cnt_pre_event", set_cycle_db_count);
    reg("cycle_db_positions_pre_event", set_cycle_db_positions);
    reg("cycle_db_positions_post_event", check_cycle_db_positions);
    reg("cycle_db_position_number_pre_event", set_cycle_db_position);
    reg("cycle_db_position_cis_trans_pre_event", set_cycle_db_position_cistrans);
    reg("cylce_element_pre_event", add_cycle_element);
    
    // set linkage events
    reg("fatty_acyl_linkage_pre_event", set_acyl_linkage);
    reg("fatty_acyl_linkage_post_event", add_acyl_linkage);
    reg("fatty_alkyl_linkage_pre_event", set_alkyl_linkage);
    reg("fatty_alkyl_linkage_post_event", add_alkyl_linkage);
    reg("fatty_linkage_number_pre_event", set_fatty_linkage_number);
    reg("fatty_acyl_linkage_sign_pre_event", set_linkage_type);
    reg("hydrocarbon_chain_pre_event", set_hydrocarbon_chain);
    reg("hydrocarbon_chain_post_event", add_hydrocarbon_chain);
    reg("hydrocarbon_number_pre_event", set_fatty_linkage_number);
    
    // set remaining events
    reg("ring_stereo_pre_event", set_ring_stereo);
    reg("pl_hg_fa_pre_event", set_hg_acyl);
    reg("pl_hg_fa_post_event", add_hg_acyl);
    reg("pl_hg_alk_pre_event", set_hg_alkyl);
    reg("pl_hg_alk_post_event", add_hg_alkyl);
    reg("pl_hg_species_pre_event", add_pl_species_data);
    reg("hg_pip_m_pre_event", suffix_decorator_molecular);
    reg("hg_pip_d_pre_event", suffix_decorator_molecular);
    reg("hg_pip_t_pre_event", suffix_decorator_molecular);
    reg("hg_PE_PS_type_pre_event", suffix_decorator_species);
    reg("acer_hg_post_event", set_acer);
    reg("acer_species_post_event", set_acer_species);
    
    reg("sterol_definition_post_event", set_sterol_definition);
    reg("adduct_heavy_element_pre_event", set_heavy_element);
    reg("adduct_heavy_number_pre_event", set_heavy_number);
    reg("adduct_heavy_component_post_event", add_heavy_component);
    
    debug = "";
}


ShorthandParserEventHandler::~ShorthandParserEventHandler(){
}


const set<string> ShorthandParserEventHandler::special_types {"acyl", "alkyl", "decorator_acyl", "decorator_alkyl", "cc"};



void ShorthandParserEventHandler::reset_lipid(TreeNode *node) {
    level = COMPLETE_STRUCTURE;
    adduct = NULL;
    head_group = "";
    fa_list->clear();
    current_fas.clear();
    headgroup_decorators->clear();
    tmp.remove_all();
    acer_species = false;
    contains_stereo_information = false;
    heavy_element = ELEMENT_C;
    heavy_element_number = 0;
}


void ShorthandParserEventHandler::set_sterol_definition(TreeNode *node){
    head_group += " " + node->get_text();
    fa_list->erase(fa_list->begin());
}



void ShorthandParserEventHandler::set_carbohydrate_number(TreeNode *node){
    int carbohydrate_num = node->get_int();
    if (!headgroup_decorators->empty() && carbohydrate_num > 0){
        headgroup_decorators->back()->count += (carbohydrate_num - 1);
    }
}
            
            

void ShorthandParserEventHandler::set_glyco_sphingo_lipid(TreeNode *node){
    
}


void ShorthandParserEventHandler::build_lipid(TreeNode *node) {
    if (acer_species) fa_list->at(0)->num_carbon -= 2;
    Headgroup *headgroup = prepare_headgroup_and_checks(false);
    
    // add count numbers for fatty acyl chains
    int fa_it = !fa_list->empty() && (fa_list->front()->lipid_FA_bond_type == LCB_REGULAR || fa_list->front()->lipid_FA_bond_type == LCB_EXCEPTION);
    for (int it = fa_it; it < (int)fa_list->size(); ++it){
        fa_list->at(it)->name += std::to_string(it + 1);
    }
    
    LipidAdduct *lipid = new LipidAdduct();
    lipid->adduct = adduct;
    lipid->lipid = assemble_lipid(headgroup);
    
    if (tmp.contains_key("num_ethers")) lipid->lipid->info->num_ethers = tmp.get_int("num_ethers");
    
    BaseParserEventHandler<LipidAdduct*>::content = lipid;
}


void ShorthandParserEventHandler::set_acer(TreeNode *node){
    head_group = "ACer";
    HeadgroupDecorator *hgd = new HeadgroupDecorator("decorator_acyl", -1, 1, 0, true);
    hgd->functional_groups->insert({"decorator_acyl", vector<FunctionalGroup*>{fa_list->back()}});
    fa_list->pop_back();
    headgroup_decorators->push_back(hgd);
}
    
    
void ShorthandParserEventHandler::set_acer_species(TreeNode *node){
    head_group = "ACer";
    set_lipid_level(SPECIES);
    HeadgroupDecorator *hgd = new HeadgroupDecorator("decorator_acyl", -1, 1, 0, true);
    hgd->functional_groups->insert({"decorator_acyl", vector<FunctionalGroup*>{new FattyAcid("FA", 2)}});
    headgroup_decorators->push_back(hgd);
    acer_species = true;
}



void ShorthandParserEventHandler::add_cycle_element(TreeNode *node){
    string element = node->get_text();
    
    if (uncontains_val(element_positions, element)){
        throw LipidParsingException("Element '" + element + "' unknown");
    }
    
    tmp.get_dictionary(FA_I)->get_list("cycle_elements")->add_int(element_positions.at(element));
}



void ShorthandParserEventHandler::set_headgroup_name(TreeNode *node){
    if (head_group.size() == 0) head_group = node->get_text();
}



void ShorthandParserEventHandler::set_carbohydrate(TreeNode *node){
    string carbohydrate = node->get_text();
    FunctionalGroup* functional_group = 0;
    try {
        functional_group = KnownFunctionalGroups::get_functional_group(carbohydrate);
    }
    catch (const std::exception& e){
        throw LipidParsingException("Carbohydrate '" + carbohydrate + "' unknown");
    }
    
    functional_group->elements->at(ELEMENT_O) -= 1;
    if (tmp.contains_key("func_group_head") && tmp.get_int("func_group_head") == 1){
        headgroup_decorators->push_back((HeadgroupDecorator*)functional_group);
    }
    else {
        if (uncontains_val_p(current_fas.at(current_fas.size() - 1)->functional_groups, carbohydrate)){
            current_fas.at(current_fas.size() - 1)->functional_groups->insert({carbohydrate, vector<FunctionalGroup*>()});
        }
        current_fas.at(current_fas.size() - 1)->functional_groups->at(carbohydrate).push_back(functional_group);
    }
}



void ShorthandParserEventHandler::set_carbohydrate_sn_position(TreeNode *node){
    set_lipid_level(SN_POSITION);
    tmp.set_int("func_group_head", 1);
}



void ShorthandParserEventHandler::set_carbohydrate_isomeric(TreeNode *node){
    tmp.set_int("func_group_head", 1);
}



void ShorthandParserEventHandler::suffix_decorator_molecular(TreeNode *node){
    headgroup_decorators->push_back(new HeadgroupDecorator(node->get_text(), -1, 1, 0, true, MOLECULAR_SPECIES));
}



void ShorthandParserEventHandler::suffix_decorator_species(TreeNode *node){
    headgroup_decorators->push_back(new HeadgroupDecorator(node->get_text(), -1, 1, 0, true, SPECIES));
}



void ShorthandParserEventHandler::set_pl_hg_triple(TreeNode *node){
        set_molecular_level(node);
        set_headgroup_name(node);
}



void ShorthandParserEventHandler::set_ring_stereo(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_ring_stereo", node->get_text());
}



void ShorthandParserEventHandler::set_fatty_acyl_stereo(TreeNode *node){
    current_fas.back()->stereochemistry = node->get_text();
    contains_stereo_information = true;
}


void ShorthandParserEventHandler::add_dihydroxyl(TreeNode* node){
    if (uncontains_val(LCB_STATES, ((FattyAcid*)current_fas.back())->lipid_FA_bond_type)) return;
    
    int num_h = 1;
    if (contains_val(SP_EXCEPTION_CLASSES, head_group) && headgroup_decorators->size() == 0) num_h += 1;
    
    FunctionalGroup* functional_group = KnownFunctionalGroups::get_functional_group("OH");
    functional_group->count = num_h;
    if (uncontains_val_p(current_fas.back()->functional_groups, "OH")) current_fas.back()->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
    current_fas.back()->functional_groups->at("OH").push_back(functional_group);
}



void ShorthandParserEventHandler::add_pl_species_data(TreeNode *node){
    set_lipid_level(SPECIES);
    HeadgroupDecorator *hgd = new HeadgroupDecorator("");
    hgd->elements->at(ELEMENT_O) += 1;
    hgd->elements->at(ELEMENT_H) -= 1;
    headgroup_decorators->push_back(hgd);
}



void ShorthandParserEventHandler::new_fatty_acyl_chain(TreeNode *node){
    current_fas.push_back(new FattyAcid("FA"));
    tmp.set_dictionary(FA_I, new GenericDictionary());
}



void ShorthandParserEventHandler::new_lcb(TreeNode *node){
    new_fatty_acyl_chain(node);
    ((FattyAcid*)current_fas.back())->set_type(LCB_REGULAR);
    current_fas.back()->name = "LCB";
}



void ShorthandParserEventHandler::add_fatty_acyl_chain(TreeNode *node){
    string fg_i = "fa" + std::to_string(current_fas.size() - 2);
    string special_type = "";
    if (current_fas.size() >= 2 && tmp.contains_key(fg_i) && tmp.get_dictionary(fg_i)->contains_key("fg_name")){
        string fg_name = tmp.get_dictionary(fg_i)->get_string("fg_name");
        if (contains_val(special_types, fg_name)){
            special_type = fg_name;
        }
    }
    
    string fa_i = FA_I;
    if (current_fas.back()->double_bonds->get_num() != tmp.get_dictionary(fa_i)->get_int("db_count")){
        throw LipidException("Double bond count does not match with number of double bond positions");
    }
    else if (current_fas.back()->double_bonds->get_num() > 0 && current_fas.back()->double_bonds->double_bond_positions.size() == 0){
        set_lipid_level(SN_POSITION);
    }
    tmp.remove(fa_i);
    
    FattyAcid* fa = (FattyAcid*)current_fas.back();
    current_fas.pop_back();
    if (special_type.length() > 0){
        fa->name = special_type;
        if (uncontains_val_p(current_fas.back()->functional_groups, special_type)) current_fas.back()->functional_groups->insert({special_type, vector<FunctionalGroup*>()});
        current_fas.back()->functional_groups->at(special_type).push_back(fa);
    }
    else {
        fa_list->push_back(fa);
    }
}



void ShorthandParserEventHandler::set_double_bond_position(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_int("db_position", node->get_int());
}



void ShorthandParserEventHandler::set_double_bond_information(TreeNode *node){
    string fa_i = FA_I;
    tmp.get_dictionary(fa_i)->set_int("db_position", 0);
    tmp.get_dictionary(fa_i)->set_string("db_cistrans", "");
}



void ShorthandParserEventHandler::add_double_bond_information(TreeNode *node){
    string fa_i = FA_I;
    int pos = tmp.get_dictionary(fa_i)->get_int("db_position");
    string cistrans = tmp.get_dictionary(fa_i)->get_string("db_cistrans");
    
    if (cistrans == ""){
        set_lipid_level(STRUCTURE_DEFINED);
    }
    
    tmp.get_dictionary(fa_i)->remove("db_position");
    tmp.get_dictionary(fa_i)->remove("db_cistrans");
    current_fas.back()->double_bonds->double_bond_positions.insert({pos, cistrans});
}



void ShorthandParserEventHandler::set_cistrans(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("db_cistrans", node->get_text());
}



void ShorthandParserEventHandler::set_functional_group(TreeNode *node){
    string fa_i = FA_I;
    GenericDictionary* gd = tmp.get_dictionary(fa_i);
    gd->set_int("fg_pos", -1);
    gd->set_string("fg_name", "0");
    gd->set_int("fg_cnt", 1);
    gd->set_string("fg_stereo", "");
    gd->set_string("fg_ring_stereo", "");
}



void ShorthandParserEventHandler::set_cycle(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_name", "cy");
    current_fas.push_back(new Cycle(0));
    
    string fa_i = FA_I;
    tmp.set_dictionary(fa_i, new GenericDictionary());
    tmp.get_dictionary(fa_i)->set_list("cycle_elements", new GenericList());
}



void ShorthandParserEventHandler::add_cycle(TreeNode *node){
    string fa_i = FA_I;
    GenericList *cycle_elements = tmp.get_dictionary(fa_i)->get_list("cycle_elements");
    Cycle *cycle = (Cycle*)current_fas.back();
    current_fas.pop_back();
    for (int i = 0; i < (int)cycle_elements->list.size(); ++i){
        cycle->bridge_chain->push_back((Element)cycle_elements->get_int(i));
    }
    tmp.get_dictionary(fa_i)->remove("cycle_elements");
        
    if (cycle->start > -1 && cycle->end > -1 && cycle->end - cycle->start + 1 + (int)cycle->bridge_chain->size() < cycle->cycle){
        throw ConstraintViolationException("Cycle length '" + std::to_string(cycle->cycle) + "' does not match with cycle description.");
    }
    if (uncontains_val_p(current_fas.back()->functional_groups, "cy")){
        current_fas.back()->functional_groups->insert({"cy", vector<FunctionalGroup*>()});
    }
    current_fas.back()->functional_groups->at("cy").push_back(cycle);
}



void ShorthandParserEventHandler::set_fatty_linkage_number(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_int("linkage_pos", node->get_int());
}



void ShorthandParserEventHandler::set_hg_acyl(TreeNode *node){
    string fa_i = FA_I;
    tmp.set_dictionary(fa_i, new GenericDictionary());
    tmp.get_dictionary(fa_i)->set_string("fg_name", "decorator_acyl");
    current_fas.push_back(new HeadgroupDecorator("decorator_acyl", -1, 1, 0, true));
    tmp.set_dictionary(FA_I, new GenericDictionary());
}



void ShorthandParserEventHandler::add_hg_acyl(TreeNode *node){
    tmp.remove(FA_I);
    headgroup_decorators->push_back((HeadgroupDecorator*)current_fas.back());
    current_fas.pop_back();
    tmp.remove(FA_I);
}



void ShorthandParserEventHandler::set_hg_alkyl(TreeNode *node){
    tmp.set_dictionary(FA_I, new GenericDictionary());
    tmp.get_dictionary(FA_I)->set_string("fg_name", "decorator_alkyl");
    current_fas.push_back(new HeadgroupDecorator("decorator_alkyl", -1, 1, 0, true));
    tmp.set_dictionary(FA_I, new GenericDictionary());
}



void ShorthandParserEventHandler::add_hg_alkyl(TreeNode *node){
    tmp.remove(FA_I);
    headgroup_decorators->push_back((HeadgroupDecorator*)current_fas.back());
    current_fas.pop_back();
    tmp.remove(FA_I);
}



void ShorthandParserEventHandler::set_linkage_type(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_int("linkage_type", node->get_text() == "N");
}



void ShorthandParserEventHandler::set_hydrocarbon_chain(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_name", "cc");
    current_fas.push_back(new CarbonChain((FattyAcid*)0));
    tmp.set_dictionary(FA_I, new GenericDictionary());
    tmp.get_dictionary(FA_I)->set_int("linkage_pos", -1);
}



void ShorthandParserEventHandler::add_hydrocarbon_chain(TreeNode *node){
    int linkage_pos = tmp.get_dictionary(FA_I)->get_int("linkage_pos");
    tmp.remove(FA_I);
    CarbonChain *cc = (CarbonChain*)current_fas.back();
    current_fas.pop_back();
    cc->position = linkage_pos;
    if (linkage_pos == -1) set_lipid_level(STRUCTURE_DEFINED);
    
    if (uncontains_val_p(current_fas.back()->functional_groups, "cc")) current_fas.back()->functional_groups->insert({"cc", vector<FunctionalGroup*>()});
    current_fas.back()->functional_groups->at("cc").push_back(cc);
}



void ShorthandParserEventHandler::set_acyl_linkage(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_name", "acyl");
    current_fas.push_back(new AcylAlkylGroup((FattyAcid*)0));
    tmp.set_dictionary(FA_I, new GenericDictionary());
    tmp.get_dictionary(FA_I)->set_int("linkage_pos", -1);
}



void ShorthandParserEventHandler::add_acyl_linkage(TreeNode *node){
    bool linkage_type = tmp.get_dictionary(FA_I)->get_int("linkage_type");
    int linkage_pos = tmp.get_dictionary(FA_I)->get_int("linkage_pos");
    
    tmp.remove(FA_I);
    AcylAlkylGroup *acyl = (AcylAlkylGroup*)current_fas.back();
    current_fas.pop_back();
        
    acyl->position = linkage_pos;
    acyl->set_N_bond_type(linkage_type);
    if (linkage_pos == -1) set_lipid_level(STRUCTURE_DEFINED);
        
    if (uncontains_val_p(current_fas.back()->functional_groups, "acyl")) current_fas.back()->functional_groups->insert({"acyl", vector<FunctionalGroup*>()});
    current_fas.back()->functional_groups->at("acyl").push_back(acyl);
}



void ShorthandParserEventHandler::set_alkyl_linkage(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_name", "alkyl");
    current_fas.push_back(new AcylAlkylGroup(0, -1, 1, true));
    tmp.set_dictionary(FA_I, new GenericDictionary());
    tmp.get_dictionary(FA_I)->set_int("linkage_pos", -1);
}



void ShorthandParserEventHandler::add_alkyl_linkage(TreeNode *node){
    int linkage_pos = tmp.get_dictionary(FA_I)->get_int("linkage_pos");
    tmp.remove(FA_I);
    AcylAlkylGroup *alkyl = (AcylAlkylGroup*)current_fas.back();
    current_fas.pop_back();
    
    alkyl->position = linkage_pos;
    if (linkage_pos == -1) set_lipid_level(STRUCTURE_DEFINED);
    
    if (uncontains_val_p(current_fas.back()->functional_groups, "alkyl")) current_fas.back()->functional_groups->insert({"alkyl", vector<FunctionalGroup*>()});
    current_fas.back()->functional_groups->at("alkyl").push_back(alkyl);
}



void ShorthandParserEventHandler::set_cycle_start(TreeNode *node){
    ((Cycle*)current_fas.back())->start = node->get_int();
    ((Cycle*)current_fas.back())->position = node->get_int();
}



void ShorthandParserEventHandler::set_cycle_end(TreeNode *node){
    ((Cycle*)current_fas.back())->end = node->get_int();
}



void ShorthandParserEventHandler::set_cycle_number(TreeNode *node){
    ((Cycle*)current_fas.back())->cycle = node->get_int();
}



void ShorthandParserEventHandler::set_cycle_db_count(TreeNode *node){
    ((Cycle*)current_fas.back())->double_bonds->num_double_bonds = node->get_int();
}



void ShorthandParserEventHandler::set_cycle_db_positions(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_int("cycle_db", ((Cycle*)current_fas.back())->double_bonds->get_num());
}



void ShorthandParserEventHandler::check_cycle_db_positions(TreeNode *node){
    if (((Cycle*)current_fas.back())->double_bonds->get_num() != tmp.get_dictionary(FA_I)->get_int("cycle_db")){
        throw LipidException("Double bond number in cycle does not correspond to number of double bond positions.");
    }
}



void ShorthandParserEventHandler::set_cycle_db_position(TreeNode *node){
    int pos = node->get_int();
    ((Cycle*)current_fas.back())->double_bonds->double_bond_positions.insert({pos, ""});
    tmp.get_dictionary(FA_I)->set_int("last_db_pos", pos);
}



void ShorthandParserEventHandler::set_cycle_db_position_cistrans(TreeNode *node){
    int pos = tmp.get_dictionary(FA_I)->get_int("last_db_pos");
    ((Cycle*)current_fas.back())->double_bonds->double_bond_positions.at(pos) = node->get_text();
}



void ShorthandParserEventHandler::set_functional_group_position(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_int("fg_pos", node->get_int());
}



void ShorthandParserEventHandler::set_functional_group_name(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_name", node->get_text());
}



void ShorthandParserEventHandler::set_functional_group_count(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_int("fg_cnt", node->get_int());
}



void ShorthandParserEventHandler::set_functional_group_stereo(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_stereo", node->get_text());
    contains_stereo_information = true;
}



void ShorthandParserEventHandler::set_sn_position_func_group(TreeNode *node){
    tmp.get_dictionary(FA_I)->set_string("fg_name", node->get_text());
    set_lipid_level(SN_POSITION);
}


void ShorthandParserEventHandler::add_functional_group(TreeNode *node){
    string fa_i = FA_I;
    GenericDictionary *gd = tmp.get_dictionary(FA_I);
    string fg_name = gd->get_string("fg_name");
    
    if (contains_val(special_types, fg_name) || fg_name == "cy") return;
        
    int fg_pos = gd->get_int("fg_pos");
    int fg_cnt = gd->get_int("fg_cnt");
    string fg_stereo = gd->get_string("fg_stereo");
    string fg_ring_stereo = gd->get_string("fg_ring_stereo");
    
    if (fg_pos == -1){
        set_lipid_level(STRUCTURE_DEFINED);
    }
    
    if (fg_cnt <= 0){
        return;
    }
    
    FunctionalGroup *functional_group = 0;
    try {
        functional_group = KnownFunctionalGroups::get_functional_group(fg_name);
    }
    catch (const std::exception& e) {
        throw LipidParsingException("'" + fg_name + "' unknown");
    }
    
    functional_group->position = fg_pos;
    functional_group->count = fg_cnt;
    functional_group->stereochemistry = fg_stereo;
    functional_group->ring_stereo = fg_ring_stereo;
    
    
    gd->remove("fg_pos");
    gd->remove("fg_name");
    gd->remove("fg_cnt");
    gd->remove("fg_stereo");
    
    if (uncontains_val_p(current_fas.back()->functional_groups, fg_name)) current_fas.back()->functional_groups->insert({fg_name, vector<FunctionalGroup*>()});
    current_fas.back()->functional_groups->at(fg_name).push_back(functional_group);
}



void ShorthandParserEventHandler::set_ether_type(TreeNode *node){
    string ether_type = node->get_text();
    if (ether_type == "O-") ((FattyAcid*)current_fas.back())->lipid_FA_bond_type = ETHER_PLASMANYL;
    else if (ether_type == "P-") ((FattyAcid*)current_fas.back())->lipid_FA_bond_type = ETHER_PLASMENYL;
}



void ShorthandParserEventHandler::set_ether_num(TreeNode *node){
    int num_ethers = 0;
    string ether = node->get_text();
    if (ether == "d") num_ethers = 2;
    else if (ether == "t") num_ethers = 3;
    else if (ether == "e") num_ethers = 4;
    tmp.set_int("num_ethers", num_ethers);
}



void ShorthandParserEventHandler::set_species_level(TreeNode *node){
    set_lipid_level(SPECIES);
}



void ShorthandParserEventHandler::set_molecular_level(TreeNode *node){
    set_lipid_level(MOLECULAR_SPECIES);
}



void ShorthandParserEventHandler::set_carbon(TreeNode *node){
    ((FattyAcid*)current_fas.back())->num_carbon = node->get_int();
}



void ShorthandParserEventHandler::set_double_bond_count(TreeNode *node){
    int db_cnt = node->get_int();
    tmp.get_dictionary(FA_I)->set_int("db_count", db_cnt);
    ((FattyAcid*)current_fas.back())->double_bonds->num_double_bonds = db_cnt;
}



void ShorthandParserEventHandler::new_adduct(TreeNode *node){
    if (!adduct) adduct = new Adduct("", "");
}



void ShorthandParserEventHandler::add_adduct(TreeNode *node){
    adduct->adduct_string = node->get_text();
}



void ShorthandParserEventHandler::add_charge(TreeNode *node){
    adduct->charge = node->get_int();
}



void ShorthandParserEventHandler::add_charge_sign(TreeNode *node){
    string sign = node->get_text();
    if (sign == "+") adduct->set_charge_sign(1);
    else adduct->set_charge_sign(-1);
    if (adduct->charge == 0) adduct->charge = 1;
}



void ShorthandParserEventHandler::set_heavy_element(TreeNode *node){
    heavy_element = heavy_element_table.at(node->get_text());
}



void ShorthandParserEventHandler::set_heavy_number(TreeNode *node){
    heavy_element_number = node->get_int();
}



void ShorthandParserEventHandler::add_heavy_component(TreeNode *node){
    adduct->heavy_elements[heavy_element] = heavy_element_number;
}



    
