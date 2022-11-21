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

const map<string, vector<string> > LipidBaseParserEventHandler::glyco_table{{"ga1", {"Gal", "GalNAc", "Gal", "Glc"}},
               {"ga2", {"GalNAc", "Gal", "Glc"}},
               {"gb3", {"Gal", "Gal", "Glc"}},
               {"gb4", {"GalNAc", "Gal", "Gal", "Glc"}},
               {"gd1", {"Gal", "GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gd1a", {"Hex", "Hex", "Hex", "HexNAc", "NeuAc", "NeuAc"}},
               {"gd2", {"GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gd3", {"NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gm1", {"Gal", "GalNAc", "NeuAc", "Gal", "Glc"}},
               {"gm2", {"GalNAc", "NeuAc", "Gal", "Glc"}},
               {"gm3", {"NeuAc", "Gal", "Glc"}},
               {"gm4", {"NeuAc", "Gal"}},
               {"gp1", {"NeuAc", "NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gq1", {"NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gt1", {"Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gt2", {"GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
               {"gt3", {"NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}}};


LipidBaseParserEventHandler::LipidBaseParserEventHandler() : BaseParserEventHandler<LipidAdduct*>() {
    fa_list = new vector<FattyAcid*>();
    level = FULL_STRUCTURE;
    head_group = "";
    lcb = NULL;
    current_fa = NULL;
    adduct = NULL;
    headgroup_decorators = new vector<HeadgroupDecorator*>();
    adduct = NULL;
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



Headgroup* LipidBaseParserEventHandler::prepare_headgroup_and_checks(){
    
    string hg = to_lower(head_group);
    if (contains_val(glyco_table, hg)){
    
        for (auto carbohydrate : glyco_table.at(hg)){
            FunctionalGroup* functional_group = 0;
            try {
                functional_group = KnownFunctionalGroups::get_functional_group(carbohydrate);
            }
            catch (const std::exception& e){
                throw LipidParsingException("Carbohydrate '" + carbohydrate + "' unknown");
            }
            
            functional_group->elements->at(ELEMENT_O) -= 1;
            headgroup_decorators->push_back((HeadgroupDecorator*)functional_group);
        }
        head_group = "Cer";
    }
    
    
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
        
