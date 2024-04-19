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


#include "cppgoslin/domain/LipidSpeciesInfo.h"

LipidSpeciesInfo::LipidSpeciesInfo (LipidClass _lipid_class) : FattyAcid("info") {
    lipid_class = _lipid_class;
    level = NO_LEVEL;
    num_ethers = 0;
    num_specified_fa = 0;
    extended_class = ESTER;
    ClassMap &lipid_classes = LipidClasses::get_instance().lipid_classes;
    total_fa = contains_val(lipid_classes, lipid_class) ? lipid_classes.at(lipid_class).max_num_fa : 0;
    poss_fa = contains_val(lipid_classes, lipid_class) ? lipid_classes.at(lipid_class).possible_num_fa : 0;
}



LipidSpeciesInfo *LipidSpeciesInfo::copy(){
    LipidSpeciesInfo *lsi = new LipidSpeciesInfo(lipid_class);
    lsi->level = level;
    lsi->num_ethers = num_ethers;
    lsi->num_specified_fa = num_specified_fa;
    lsi->position = position;
    lsi->total_fa = total_fa;
    lsi->extended_class = extended_class;
    lsi->num_carbon = num_carbon;
    delete lsi->double_bonds;
    lsi->double_bonds = double_bonds->copy();
    lsi->lipid_FA_bond_type = lipid_FA_bond_type;
    
    for (auto kv : *functional_groups){
        lsi->functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
        for (auto func_group : kv.second){
            lsi->functional_groups->at(kv.first).push_back(func_group->copy());
        }
    }
    return lsi;
}



ElementTable* LipidSpeciesInfo::get_elements(){
    ElementTable* elements = FattyAcid::get_elements();
    if (lipid_FA_bond_type != LCB_EXCEPTION) elements->at(ELEMENT_O) -= (num_ethers == 0);
    elements->at(ELEMENT_H) += num_ethers == 0 ? 1 : -1;
    
    return elements;
}


void LipidSpeciesInfo::add(FattyAcid* _fa){
    if ((_fa->lipid_FA_bond_type == ETHER_PLASMENYL || _fa->lipid_FA_bond_type == ETHER_PLASMANYL) && _fa->lipid_FA_bond_type != LCB_EXCEPTION && _fa->lipid_FA_bond_type != LCB_REGULAR){
        num_ethers += 1;
        lipid_FA_bond_type = ETHER_PLASMANYL;
        extended_class = _fa->lipid_FA_bond_type;
    }
    else if (_fa->lipid_FA_bond_type == LCB_EXCEPTION || _fa->lipid_FA_bond_type == LCB_REGULAR){
        lipid_FA_bond_type = _fa->lipid_FA_bond_type;
    }
            
    else{
        num_specified_fa += 1;
    }
    for (auto &kv : *(_fa->functional_groups)){
        if (uncontains_val_p(functional_groups, kv.first)) functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
        for (auto func_group : kv.second) {
            functional_groups->at(kv.first).push_back(func_group->copy());
        }
    }
        
    ElementTable* e = _fa->get_elements();
    num_carbon += e->at(ELEMENT_C);
    delete e;
    double_bonds->num_double_bonds += _fa->get_double_bonds();
}


string LipidSpeciesInfo::to_string(){
    stringstream info_string;
    info_string << ether_prefix[num_ethers];
    info_string << num_carbon << ":" << double_bonds->get_num();
    
    
    ElementTable *elements = get_functional_group_elements();
    for (int i = 2; i < (int)element_order.size(); ++i){
        Element e = element_order.at(i);
        if (elements->at(e) > 0){
            info_string << ";" << element_shortcut.at(e);
            if (elements->at(e) > 1){
                info_string << elements->at(e);
            }
        }
    } 
    delete elements;
    
        
    return info_string.str();
}
