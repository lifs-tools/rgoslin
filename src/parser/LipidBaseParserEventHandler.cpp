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


#include "cppgoslin/parser/LipidBaseParserEventHandler.h"
               
const map<string, int> LipidBaseParserEventHandler::fa_synonyms{{"Palmitic acid", 0}, {"Linoleic acid", 1}, {"AA", 2}, {"ALA", 3}, {"EPA", 4}, {"DHA", 5}, {"LTB4", 6}, {"Resolvin D3", 7}, {"Maresin 1", 8},  {"Resolvin D2", 9}, {"Resolvin D5", 10}, {"Resolvin D1", 11}, {"TXB1", 12}, {"TXB2", 13}, {"TXB3", 14}, {"PGF2alpha", 15}, {"PGD2", 16}, {"PGE2", 17}, {"PGB2", 18}, {"15d-PGJ2", 19}, {"PGJ2", 20}, {"Arachidonic acid", 21}, {"LA", 22}, {"Mar1", 23}, {"RvD3", 24}, {"PGF1alpha", 25}, {"PDX", 26}, {"Oleic acid", 27}, {"OA", 28}, {"DGLA", 29}, {"iPF2alpha-VI", 30}};

LipidBaseParserEventHandler::LipidBaseParserEventHandler() : BaseParserEventHandler<LipidAdduct*>() {
    fa_list = new vector<FattyAcid*>();
    level = FULL_STRUCTURE;
    head_group = "";
    lcb = NULL;
    current_fa = NULL;
    adduct = NULL;
    headgroup_decorators = new vector<HeadgroupDecorator*>();
    use_head_group = false;
}

const set<string> LipidBaseParserEventHandler::SP_EXCEPTION_CLASSES{"Cer", "Ceramide", "Sphingosine", "So", "Sphinganine", "Sa", "SPH", "Sph", "LCB"};



LipidBaseParserEventHandler::~LipidBaseParserEventHandler(){
    delete fa_list;
    delete headgroup_decorators;
}



void LipidBaseParserEventHandler::set_lipid_level(LipidLevel _level){
    level = min(level, _level);
}



bool LipidBaseParserEventHandler::sp_regular_lcb(){
    return Headgroup::get_category(head_group) == SP && contains_val(LCB_STATES, current_fa->lipid_FA_bond_type) && !(contains_val(LipidBaseParserEventHandler::SP_EXCEPTION_CLASSES, head_group) && headgroup_decorators->size() == 0);
    
}

bool LipidBaseParserEventHandler::check_full_structure(FunctionalGroup *obj){
    bool full = true;
    
    auto* is_fa = dynamic_cast<FattyAcid*>(obj);
    if (is_fa && ((FattyAcid*)obj)->num_carbon == 0) return true;
    if (is_fa && obj->double_bonds->num_double_bonds > 0 && obj->double_bonds->double_bond_positions.empty()) return false;
    if (is_fa && !obj->double_bonds->double_bond_positions.empty()){
        int sum = 0;
        for (auto &kv : obj->double_bonds->double_bond_positions) sum += kv.second == "E" || kv.second == "Z" || (kv.second == "" && kv.first == ((FattyAcid*)obj)->num_carbon - 1);
        full &= sum == obj->double_bonds->double_bond_positions.size();
        
    }
    
    for (auto &kv : *(obj->functional_groups)){
        for (auto fg : kv.second){
            if (fg->name == "X") continue;
            if (fg->position < 0){
                return false;
            }
            full &= check_full_structure(fg);
        }
    }
    return full;
}


FattyAcid* LipidBaseParserEventHandler::resolve_fa_synonym(string mediator_name){
    
    if (uncontains_val(fa_synonyms, mediator_name)) return 0;
    
    switch(LipidBaseParserEventHandler::fa_synonyms.at(mediator_name)){
        case 0: // Palmitic acid
            return new FattyAcid("FA", 16);
            break;
            
        case 1: // Linoleic acid
        case 22: //LA
            return new FattyAcid("FA", 18, new DoubleBonds(2));
            break;
            
        case 2: // AA":
        case 21: // Arachidonic acid
            return new FattyAcid("FA", 20, new DoubleBonds(4));
            break;
            
        case 3: // ALA":
            return new FattyAcid("FA", 18, new DoubleBonds(3));
            break;
            
        case 4: // EPA":
            return new FattyAcid("FA", 20, new DoubleBonds(5));
            break;
            
        case 5: // DHA":
            return new FattyAcid("FA", 22, new DoubleBonds(6));
            break;
            
        case 6: // LTB4
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 5;
                f2->position = 12;
                return new FattyAcid("FA", 20, new DoubleBonds({{6, "Z"}, {8, "E"}, {10, "E"}, {14, "Z"}}), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
            }
            break;
            
        case 7: // Resolvin D3
        case 24: // RvD3
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 4;
                f2->position = 11;
                f3->position = 17;
                return new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2, f3}}});
            }
            break;
            
        case 8: // Maresin 1
        case 23: // Mar1
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 4;
                f2->position = 14;
                return new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
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
                return new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2, f3}}});
            }
            break;
            
        case 10: // Resolvin D5
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 7;
                f2->position = 17;
                return new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
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
                return new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2, f3}}});
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
                return new FattyAcid("FA", 20, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
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
                return new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
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
                return new FattyAcid("FA", 20, new DoubleBonds(3), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
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
                return new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
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
                return new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
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
                return new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 18: // PGB2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 15;
                f2->position = 9;
                Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"OH", {f2}}});
                return new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 19: // 15d-PGJ2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
                f1->position = 15;
                f2->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"oxo", {f2}}});
                return new FattyAcid("FA", 20, new DoubleBonds(3), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 20: // PGJ2
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("oxo");
                f1->position = 15;
                f2->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"oxo", {f2}}});
                return new FattyAcid("FA", 20, new DoubleBonds(2), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
            
            
        case 25: // PGF1alpha
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 15;
                f2->position = 9;
                f3->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}});
                return new FattyAcid("FA", 20, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
            
        case 26: // PDX
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 10;
                f2->position = 17;
                return new FattyAcid("FA", 22, new DoubleBonds(6), new map<string, vector<FunctionalGroup*>>{{"OH", {f1, f2}}});
            }
            break;
            
        case 27: // Oleic acid
        case 28: // OA
            return new FattyAcid("FA", 18, new DoubleBonds(1));
            break;
            
        case 29: // DGLA
            return new FattyAcid("FA", 20, new DoubleBonds(3));
            break;
            
        case 30: // iPF2alpha-VI
            {
                FunctionalGroup *f1 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f2 = KnownFunctionalGroups::get_functional_group("OH");
                FunctionalGroup *f3 = KnownFunctionalGroups::get_functional_group("OH");
                f1->position = 5;
                f2->position = 9;
                f3->position = 11;
                Cycle* cy = new Cycle(5, 8, 12, 0, new map<string, vector<FunctionalGroup*>>{{"OH", {f2, f3}}});
                return new FattyAcid("FA", 20, new DoubleBonds(1), new map<string, vector<FunctionalGroup*>>{{"OH", {f1}}, {"cy", {cy}}});
            }
            break;
    }
    
    return 0;
}



Headgroup* LipidBaseParserEventHandler::prepare_headgroup_and_checks(){

    Headgroup *headgroup = new Headgroup(head_group, headgroup_decorators, use_head_group);
    
    if (use_head_group) return headgroup;
    head_group = headgroup->get_class_name();
    int true_fa = 0;
    for (auto fa : *fa_list){
        true_fa += fa->num_carbon > 0 || fa->double_bonds->get_num() > 0;
    }
    int poss_fa = contains_val(LipidClasses::get_instance().lipid_classes, headgroup->lipid_class) ? LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class).possible_num_fa : 0;
    
    // make lyso
    bool can_be_lyso = contains_val(LipidClasses::get_instance().lipid_classes, Headgroup::get_class("L" + head_group)) ? contains_val(LipidClasses::get_instance().lipid_classes.at(Headgroup::get_class("L" + head_group)).special_cases, "Lyso") : 0;
    
    
    if ((true_fa + 1 == poss_fa || true_fa + 2 == poss_fa) && level != SPECIES && headgroup->lipid_category == GP && can_be_lyso){
        if (true_fa + 1 == poss_fa) head_group = "L" + head_group;
        else head_group = "DL" + head_group;
        headgroup->decorators->clear();
        delete headgroup;
        headgroup = new Headgroup(head_group, headgroup_decorators, use_head_group);
        poss_fa = contains_val(LipidClasses::get_instance().lipid_classes, headgroup->lipid_class) ? LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class).possible_num_fa : 0;
    }
    
    else if ((true_fa + 1 == poss_fa || true_fa + 2 == poss_fa) && level != SPECIES && headgroup->lipid_category == GL && head_group == "TG"){
        if (true_fa + 1 == poss_fa) head_group = "DG";
        else head_group = "MG";
        headgroup->decorators->clear();
        delete headgroup;
        headgroup = new Headgroup(head_group, headgroup_decorators, use_head_group);
        poss_fa = contains_val(LipidClasses::get_instance().lipid_classes, headgroup->lipid_class) ? LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class).possible_num_fa : 0;
    }
    
    // check if all functional groups have a position to be full structure
    if (is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE)){
        for (auto fa : *fa_list){
            if (!check_full_structure(fa)){
                set_lipid_level(STRUCTURE_DEFINED);
                break;
            }
        }
    }
    
    
    if (level == SPECIES){
        if (true_fa == 0 && poss_fa != 0){
            string hg_name = headgroup->headgroup;
            delete headgroup;
            throw ConstraintViolationException("No fatty acyl information lipid class '" + hg_name + "' provided.");
        }
    }
        
    else if (true_fa != poss_fa && (is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE | STRUCTURE_DEFINED))){
        string hg_name = headgroup->headgroup;
        delete headgroup;
        throw ConstraintViolationException("Number of described fatty acyl chains (" + std::to_string(true_fa) + ") not allowed for lipid class '" + hg_name + "' (having " + std::to_string(poss_fa) + " fatty aycl chains).");
    }
    
    else if (contains_val(LipidClasses::get_instance().lipid_classes.at(Headgroup::get_class(head_group)).special_cases, "Lyso") && true_fa > poss_fa){
        string hg_name = headgroup->headgroup;
        delete headgroup;
        throw ConstraintViolationException("Number of described fatty acyl chains (" + std::to_string(true_fa) + ") not allowed for lipid class '" + hg_name + "' (having " + std::to_string(poss_fa) + " fatty aycl chains).");
    }
    
    
    if (contains_val(LipidClasses::get_instance().lipid_classes, headgroup->lipid_class)){
        
        if (contains_val(LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class).special_cases, "HC")){
            fa_list->front()->lipid_FA_bond_type = ETHER;
        }
        
        if (contains_val(LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class).special_cases, "Amide")){
            for (auto fatty : *fa_list) fatty->lipid_FA_bond_type = AMIDE;
        }
        
        int max_num_fa = LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class).max_num_fa;
        if (max_num_fa != (int)fa_list->size()) set_lipid_level(MOLECULAR_SPECIES);
    }
    
    // make LBC exception
    if (fa_list->size() > 0 && headgroup->sp_exception) fa_list->front()->set_type(LCB_EXCEPTION);
    return headgroup;
}
    
        
LipidSpecies* LipidBaseParserEventHandler::assemble_lipid(Headgroup *headgroup){
    
    for (auto fa : *fa_list){
        if (fa->stereo_information_missing()){
            set_lipid_level(FULL_STRUCTURE);
            break;
        }
    }
    
    LipidSpecies *ls = NULL;
    switch (level){
        case COMPLETE_STRUCTURE: ls = new LipidCompleteStructure(headgroup, fa_list); break;
        case FULL_STRUCTURE: ls = new LipidFullStructure(headgroup, fa_list); break;
        case STRUCTURE_DEFINED: ls = new LipidStructureDefined(headgroup, fa_list); break;
        case SN_POSITION: ls = new LipidSnPosition(headgroup, fa_list); break;
        case MOLECULAR_SPECIES: ls = new LipidMolecularSpecies(headgroup, fa_list); break;
        case SPECIES: ls = new LipidSpecies(headgroup, fa_list); break;
        default: break;
    }
    return ls;
}
        
