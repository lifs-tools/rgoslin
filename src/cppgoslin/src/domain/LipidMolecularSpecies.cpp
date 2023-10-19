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


#include "cppgoslin/domain/LipidMolecularSpecies.h"


LipidMolecularSpecies::LipidMolecularSpecies (Headgroup* _headgroup, vector<FattyAcid*> *_fa) : LipidSpecies(_headgroup, _fa) {
    info->level = MOLECULAR_SPECIES;
    for (auto fatty_acid : *_fa){
        if (contains_val(fa, fatty_acid->name)){
            throw ConstraintViolationException("FA names must be unique! FA with name " + fatty_acid->name + " was already added!");
        }
        fa.insert({fatty_acid->name, fatty_acid});
        fa_list.push_back(fatty_acid);
    }
           
            
    // add 0:0 dummys
    for (int i = (int)_fa->size(); i < info->total_fa; ++i){
        FattyAcid *fatty_acid = new FattyAcid("FA" + std::to_string(i + 1));
        fatty_acid->position = -1;
        fatty_acid->unresolved_hidden_fa = (1 < (int)_fa->size() && (int)_fa->size() < info->poss_fa);
        info->add(fatty_acid);
        fa.insert({fatty_acid->name, fatty_acid});
        fa_list.push_back(fatty_acid);
    }
}



string LipidMolecularSpecies::build_lipid_subspecies_name(LipidLevel level){
    if (level == NO_LEVEL) level = MOLECULAR_SPECIES;
    
    string fa_separator = (level != MOLECULAR_SPECIES || headgroup->lipid_category == SP) ? "/" : "_";
    stringstream lipid_name;
    
    lipid_name << headgroup->get_lipid_string(level);
    
    string fa_headgroup_separator = (headgroup->lipid_category != ST) ? " " : "/";
    
    switch (level){
        case COMPLETE_STRUCTURE:
        case FULL_STRUCTURE:
        case STRUCTURE_DEFINED:
        case SN_POSITION:
            if (fa_list.size() > 0){
                lipid_name << fa_headgroup_separator;
                int i = 0;
    
                for (auto fatty_acid : fa_list){
                    if (i++ > 0) lipid_name << fa_separator;
                    lipid_name << fatty_acid->to_string(level);
                }
            }
            break;
            
        default:
            bool go_on = false;
            for (auto fatty_acid : fa_list){
                if (fatty_acid->num_carbon > 0){
                    go_on = true;
                    break;
                }
            }
            
            if (go_on){
                lipid_name << fa_headgroup_separator;
                int i = 0;
                for (auto fatty_acid : fa_list){
                    if (fatty_acid->num_carbon > 0){
                        if (i++ > 0) lipid_name << fa_separator;
                        lipid_name << fatty_acid->to_string(level);
                    }
                }
            }
            break;
    }
    return lipid_name.str();
}


LipidLevel LipidMolecularSpecies::get_lipid_level(){
    return MOLECULAR_SPECIES;
}



ElementTable* LipidMolecularSpecies::get_elements(){
    ElementTable* elements = create_empty_table();
    
    ElementTable* hg_elements = headgroup->get_elements();
    for (auto &kv : *hg_elements) elements->at(kv.first) += kv.second;
    delete hg_elements;

    
    // add elements from all fatty acyl chains
    for (auto fatty_acid : fa_list){
        ElementTable* fa_elements = fatty_acid->get_elements();
        for (auto &kv : *fa_elements) elements->at(kv.first) += kv.second;
        delete fa_elements;
    }
    
    return elements;
}


void LipidMolecularSpecies::sort_fatty_acyl_chains(){
    if (info->level > MOLECULAR_SPECIES || fa_list.size() < 2) return;
    sort(fa_list.begin(), fa_list.end(), [] (FattyAcid *fa1, FattyAcid *fa2) {
        // treat empty fatty acids individually
        if (fa1 == 0 || fa1->num_carbon == 0) return false;
        if (fa2 == 0 || fa2->num_carbon == 0) return true;
        
        if (fa1->lipid_FA_bond_type != fa2->lipid_FA_bond_type) return fa1->lipid_FA_bond_type < fa2->lipid_FA_bond_type;
        if (fa1->num_carbon != fa2->num_carbon) return fa1->num_carbon < fa2->num_carbon;
        int db1 = fa1->double_bonds->get_num();
        int db2 = fa2->double_bonds->get_num();
        if (db1 != db2) return db1 < db2;
        ElementTable *e1 = fa1->get_elements();
        ElementTable *e2 = fa2->get_elements();
        double mass1 = goslin::get_mass(e1);
        double mass2 = goslin::get_mass(e2);
        delete e1;
        delete e2;
        return mass1 < mass2;
    });
}


string LipidMolecularSpecies::get_lipid_string(LipidLevel level) {
    switch (level){
        case NO_LEVEL:
        case MOLECULAR_SPECIES:
            return build_lipid_subspecies_name(MOLECULAR_SPECIES);
    
        case SPECIES:
        case CLASS:
        case CATEGORY:
            return LipidSpecies::get_lipid_string(level);
    
        default:
            throw IllegalArgumentException("LipidMolecularSpecies does not know how to create a lipid string for level " + std::to_string(level));
    }
}
