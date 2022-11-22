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


#include "cppgoslin/domain/Cycle.h"

Cycle::Cycle(int _cycle, int _start, int _end, DoubleBonds* _double_bonds, map<string, vector<FunctionalGroup*> >* _functional_groups, vector<Element>* _bridge_chain) : FunctionalGroup("cy", -1, 1, _double_bonds, false, "", false, 0, _functional_groups){
    count = 1;
    cycle = _cycle;
    position = _start;
    start = _start;
    end = _end;
    elements->at(ELEMENT_H) = -2;
    bridge_chain = (_bridge_chain == 0) ? new vector<Element>() : _bridge_chain;
}

Cycle* Cycle::copy(){
    DoubleBonds* db = double_bonds->copy();
    map<string, vector<FunctionalGroup*> >* fg = new map<string, vector<FunctionalGroup*> >();
    for (auto &kv : *functional_groups){
        fg->insert({kv.first, vector<FunctionalGroup*>()});
        for (auto &func_group : kv.second) {
            fg->at(kv.first).push_back(func_group->copy());
        }
    }
    vector<Element>* bc = new vector<Element>();
    for (auto &e : *bridge_chain) bc->push_back(e);
    
    return new Cycle(cycle, start, end, db, fg, bc);
}


Cycle::~Cycle(){
    delete bridge_chain;
}
                
                
int Cycle::get_double_bonds(){
    return FunctionalGroup::get_double_bonds() + 1;
}

void Cycle::add_position(int pos){
    start += start >= pos;
    end += end >= pos;
    FunctionalGroup::add_position(pos);
}
    
    
void Cycle::rearrange_functional_groups(FunctionalGroup* parent, int shift){
    // put everything back into parent
    for (auto &kv : double_bonds->double_bond_positions) {
        parent->double_bonds->double_bond_positions.insert({kv.first, kv.second});
    }
    delete double_bonds;
    double_bonds = new DoubleBonds();
    
    for (auto &kv : *functional_groups){
        if (!contains_val_p(parent->functional_groups, kv.first)){
            parent->functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
        }
        parent->functional_groups->at(kv.first).insert(parent->functional_groups->at(kv.first).end(), functional_groups->at(kv.first).begin(), functional_groups->at(kv.first).end());
    }
    delete functional_groups;
    functional_groups = new map<string, vector<FunctionalGroup*> >();
    
    
    // shift the cycle
    shift_positions(shift);
    
    
    // take back what's mine
    // check double bonds
    for (auto &kv : parent->double_bonds->double_bond_positions){
        if (start <= kv.first && kv.first <= end){
            double_bonds->double_bond_positions.insert({kv.first, kv.second});
        }
    }
    double_bonds->num_double_bonds = double_bonds->double_bond_positions.size();
    
    for (auto &kv : double_bonds->double_bond_positions){
        parent->double_bonds->double_bond_positions.erase(kv.first);
    }
    parent->double_bonds->num_double_bonds = parent->double_bonds->double_bond_positions.size();
    
    // check functional groups
    set<string> remove_list;
    for (auto &kv : *(parent->functional_groups)){
        vector<int> remove_item;
        
        int i = 0;
        
        for (auto func_group : kv.second){
            if (start <= func_group->position && func_group->position <= end && func_group != this){
                if (!contains_val_p(functional_groups, kv.first)){
                    functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
                }
                functional_groups->at(kv.first).push_back(func_group);
                remove_item.push_back(i);
            }
            ++i;
        }
                
        while (!remove_item.empty()){
            int pos = remove_item.back();
            remove_item.pop_back();
            kv.second.erase(kv.second.begin() + pos);
        }
        if (kv.second.empty()) remove_list.insert(kv.first);
    }
        
    for (string fg : remove_list) parent->functional_groups->erase(fg);
}    
        
void Cycle::shift_positions(int shift){
    FunctionalGroup::shift_positions(shift);
    start += shift;
    end += shift;
    DoubleBonds* db = new DoubleBonds();
    for (auto &kv : double_bonds->double_bond_positions) db->double_bond_positions.insert({kv.first + shift, kv.second});
    db->num_double_bonds = db->double_bond_positions.size();
    delete double_bonds;
    double_bonds = db;
}
    

void Cycle::compute_elements(){
    for (auto el : element_order) elements->at(el) = 0;
    elements->at(ELEMENT_H) = -2 - 2 * double_bonds->num_double_bonds;
    
    for (auto &chain_element : *bridge_chain){
        switch(chain_element){
            case ELEMENT_C:
                elements->at(ELEMENT_C) += 1;
                elements->at(ELEMENT_H) += 2;
                break;
                
            case ELEMENT_N:
                elements->at(ELEMENT_N) += 1;
                elements->at(ELEMENT_H) += 1;
                break;
                
            case ELEMENT_P:
                elements->at(ELEMENT_P) += 1;
                elements->at(ELEMENT_H) += 1;
                break;
                
            case ELEMENT_As:
                elements->at(ELEMENT_As) += 1;
                elements->at(ELEMENT_H) += 1;
                break;
                
            case ELEMENT_O:
                elements->at(ELEMENT_O) += 1;
                break;
                
            case ELEMENT_S:
                elements->at(ELEMENT_S) += 1;
                break;
                
            default:
                throw ConstraintViolationException("Element '" + element_shortcut.at(chain_element) + "' cannot be part of a cycle bridge");
        }
        
    }
        
    // add all implicit carbon chain elements
    if (start != -1 && end != -1){
        int n = max((int)(cycle - (end - start + 1 + bridge_chain->size())), 0);
        elements->at(ELEMENT_C) += n;
        elements->at(ELEMENT_H) += n << 1;
    }
}

    
string Cycle::to_string(LipidLevel level){
    stringstream cycle_string;
    cycle_string << "[";
    if (start != -1 && is_level(level, FULL_STRUCTURE | COMPLETE_STRUCTURE)){
        cycle_string << start << "-" << end;
    }
    
    if (is_level(level , FULL_STRUCTURE | COMPLETE_STRUCTURE | STRUCTURE_DEFINED) && bridge_chain->size() > 0){
        for (auto &e : *bridge_chain) cycle_string << element_shortcut.at(e);
    }
    cycle_string << "cy" << cycle;
    
    cycle_string << ":" << double_bonds->num_double_bonds;
    
    if (is_level(level , FULL_STRUCTURE | COMPLETE_STRUCTURE | STRUCTURE_DEFINED)){
        if (double_bonds->double_bond_positions.size() > 0){
            int i = 0;
            cycle_string << "(";
            
            vector<int> db_keys;
            for (auto &kv : double_bonds->double_bond_positions) db_keys.push_back(kv.first);
            sort(db_keys.begin(), db_keys.end());
            
            for (auto key : db_keys){
                string value = double_bonds->double_bond_positions.at(key);
                if (i++ > 0) cycle_string << ",";
                if (is_level(level, FULL_STRUCTURE | COMPLETE_STRUCTURE)) cycle_string << key << value;
                else  cycle_string << key;
            }
            cycle_string << ")";
        }
    }
    
    if (is_level(level, FULL_STRUCTURE | COMPLETE_STRUCTURE)){
        vector<string> fg_names;
        for (auto &kv : *functional_groups) fg_names.push_back(kv.first);
        sort(fg_names.begin(), fg_names.end(), lower_name_sort_function);
        
        for (auto &fg : fg_names){
            vector<FunctionalGroup*>& fg_list = functional_groups->at(fg);
            if (fg_list.empty()) continue;
            
            sort(fg_list.begin(), fg_list.end(), FunctionalGroup::position_sort_function);
            int i = 0;
            cycle_string << ";";
            for (auto func_group : fg_list){
                if (i++ > 0) cycle_string << ",";
                cycle_string << func_group->to_string(level);
            }
        }
    }
    
    else if (level == STRUCTURE_DEFINED){
        vector<string> fg_names;
        for (auto &kv : *functional_groups) fg_names.push_back(kv.first);
        sort(fg_names.begin(), fg_names.end(), lower_name_sort_function);
        
        for (auto &fg : fg_names){
            vector<FunctionalGroup*> &fg_list = functional_groups->at(fg);
            if (fg_list.empty()) continue;
            
            else if (fg_list.size() == 1 && fg_list.at(0)->count == 1){
                cycle_string << ";" << fg_list.front()->to_string(level);
            }
            else {
                int fg_count = 0;
                for (auto func_group : fg_list) fg_count += func_group->count;
                if (fg_count > 1){
                    cycle_string << ";(" << fg << ")" << fg_count;
                }
                else {
                    cycle_string << ";" << fg;
                }
            }
        }
    }
                
    cycle_string << "]";
    if (level == COMPLETE_STRUCTURE && stereochemistry.length() > 0) cycle_string << "[" << stereochemistry << "]";
    
    return cycle_string.str();
}
