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


#include "cppgoslin/domain/LipidSpecies.h"

using namespace std;

LipidSpecies::LipidSpecies(Headgroup* _headgroup, vector<FattyAcid*>* _fa){
    headgroup = _headgroup;
    info = new LipidSpeciesInfo(headgroup->lipid_class);
    info->level = SPECIES;
    
    // add fatty acids
    if (_fa != 0){
        int i = 0;
        bool fa_it = (_fa->size() > 0) && (_fa->at(0)->lipid_FA_bond_type == LCB_EXCEPTION || _fa->at(0)->lipid_FA_bond_type == LCB_REGULAR);
        for (auto fatty_acid : *_fa){
            fatty_acid->name = (fa_it && i == 0) ? "LCB" : "FA" + std::to_string(i + 1 - fa_it);
            fatty_acid->position = -1;
            info->add(fatty_acid);
            ++i;
        }
    }
}


LipidSpecies::~LipidSpecies(){
    for (auto fa : fa_list){
        delete fa;
    }
    
    delete info;
    delete headgroup;
}


LipidLevel LipidSpecies::get_lipid_level(){
    return SPECIES;
}


string LipidSpecies::get_lipid_string(LipidLevel level){
    switch (level){
            
        default:
            throw RuntimeException("LipidSpecies does not know how to create a lipid string for level " + std::to_string(level));
            
        case UNDEFINED_LEVEL:
            throw RuntimeException("LipidSpecies does not know how to create a lipid string for level " + std::to_string(level));
            
        case CLASS:
        case CATEGORY:
            return headgroup->get_lipid_string(level);
            
        case NO_LEVEL:
        case SPECIES:
            stringstream lipid_string;
            lipid_string << headgroup->get_lipid_string(level);
            
            if (info->elements->at(ELEMENT_C) > 0 || info->num_carbon > 0){
                
                LipidSpeciesInfo* lsi = info->copy();
                for (auto decorator : *headgroup->decorators){
                    if (decorator->name == "decorator_alkyl" || decorator->name == "decorator_acyl"){
                        ElementTable* e = decorator->get_elements();
                        lsi->num_carbon += e->at(ELEMENT_C);
                        delete e;
                        lsi->double_bonds->num_double_bonds += decorator->get_double_bonds();
                    }
                }
                
                lipid_string << (headgroup->lipid_category != ST ? " " : "/") << lsi->to_string();
                delete lsi;
            }
            return lipid_string.str();
    }
    return "";
}



void LipidSpecies::sort_fatty_acyl_chains(){
    
}



string LipidSpecies::get_extended_class(){
    LipidClassMeta &meta = LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class);
    bool special_case = (info->num_carbon > 0) ? contains_val(meta.special_cases, "Ether") : false;
    
    string class_name = info->level >= STRUCTURE_DEFINED ? headgroup->get_lipid_string(STRUCTURE_DEFINED) : headgroup->get_lipid_string(SPECIES);
    
    if (class_name == "UNDEFINED") return class_name;
    
    if (special_case){
        if(info->extended_class == ETHER_PLASMANYL || info->extended_class == ETHER_UNSPECIFIED){
            return class_name + " O";
        }
        else if (info->extended_class == ETHER_PLASMENYL){
            return class_name + " P";
        }
    }
    
    return class_name;
}


vector<FattyAcid*> LipidSpecies::get_fa_list(){
    return fa_list;
}




ElementTable* LipidSpecies::get_elements(){
    ElementTable* elements = create_empty_table();
    
    switch(info->level){
        case COMPLETE_STRUCTURE:
        case FULL_STRUCTURE:
        case STRUCTURE_DEFINED:
        case SN_POSITION:
        case MOLECULAR_SPECIES:
        case SPECIES:
            break;

        default:    
            throw LipidException("Element table cannot be computed for lipid level " + std::to_string(info->level));
    }
    if (headgroup->use_headgroup){
            throw LipidException("Element table cannot be computed for lipid level " + std::to_string(info->level));
    }
    
    
    ElementTable* hg_elements = headgroup->get_elements();
    for (auto &kv : *hg_elements) elements->at(kv.first) += kv.second;
    delete hg_elements;
    
    ElementTable* info_elements = info->get_elements();
    for (auto &kv : *info_elements) elements->at(kv.first) += kv.second;
    delete info_elements;
    
    // since only one FA info is provided, we have to treat this single information as
    // if we would have the complete information about all possible FAs in that lipid
    LipidClassMeta &meta = LipidClasses::get_instance().lipid_classes.at(headgroup->lipid_class);
    
    int additional_fa = meta.possible_num_fa;
    int remaining_H = meta.max_num_fa - additional_fa;
    int hydrochain = contains_val(meta.special_cases, "HC");
    
    elements->at(ELEMENT_O) -= -additional_fa + info->num_ethers + headgroup->sp_exception + hydrochain;
    elements->at(ELEMENT_H) += -additional_fa + remaining_H + 2 * info->num_ethers + 2 * hydrochain;
    
    if (contains_val(meta.special_cases, "Amide")){
        elements->at(ELEMENT_O) -= meta.max_num_fa;
        elements->at(ELEMENT_H) += meta.max_num_fa;
    }
    
    return elements;
}
