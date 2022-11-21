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



#include "cppgoslin/domain/FunctionalGroup.h"

FunctionalGroup::FunctionalGroup(string _name, int _position, int _count, DoubleBonds* _double_bonds, bool _is_atomic, string _stereochemistry, ElementTable* _elements, map<string, vector<FunctionalGroup*>>* _functional_groups){
    name = _name;
    position = _position;
    count = _count;
    stereochemistry = _stereochemistry;
    ring_stereo = "";
    double_bonds = (_double_bonds != 0) ? _double_bonds : new DoubleBonds(0);
    is_atomic = _is_atomic;
    num_atoms = 0;
    if (_elements != 0){
        elements = _elements;
        for (auto kv : *elements) num_atoms += kv.second;
        num_atoms = max(0, num_atoms);
    }
    else {
        elements = create_empty_table();
    }
    functional_groups = (_functional_groups != 0) ? _functional_groups : (new map<string, vector<FunctionalGroup*>>());
}


FunctionalGroup* FunctionalGroup::copy(){
    DoubleBonds* db = double_bonds->copy();
    map<string, vector<FunctionalGroup*> >* fg = new map<string, vector<FunctionalGroup*> >();
    for (auto &kv : *functional_groups){
        fg->insert({kv.first, vector<FunctionalGroup*>()});
        for (auto &func_group : kv.second) {
            fg->at(kv.first).push_back(func_group->copy());
        }
    }
    ElementTable* e = create_empty_table();
    for (auto &kv : *elements){
        e->at(kv.first) = kv.second;
    }
    
    FunctionalGroup* func_group = new FunctionalGroup(name, position, count, db, is_atomic, stereochemistry, e, fg);
    func_group->ring_stereo = ring_stereo;
    func_group->num_atoms = num_atoms;
    return func_group;
}



bool FunctionalGroup::position_sort_function (FunctionalGroup* f1, FunctionalGroup *f2) {
    return (f1->position < f2->position);
}


bool FunctionalGroup::lower_name_sort_function (string s1, string s2) {
    return (to_lower(s1) < to_lower(s2));
}


FunctionalGroup::~FunctionalGroup(){
    delete double_bonds;
    delete elements;
    for (auto &kv : *functional_groups){
        for (auto &fg : kv.second){
            delete fg;
        }
    }
    delete functional_groups;
}



ElementTable* FunctionalGroup::get_elements(){
    compute_elements();
    ElementTable* _elements = create_empty_table();
    for (auto &kv : *elements) _elements->at(kv.first) = kv.second;
    
    ElementTable* fgElements = get_functional_group_elements();
    for (auto &kv : *fgElements) _elements->at(kv.first) += kv.second;
    delete fgElements;
    return _elements;
}


void FunctionalGroup::shift_positions(int shift){
    position += shift;
    for (auto &kv : *functional_groups){
        for (auto fg : kv.second)
            fg->shift_positions(shift);
    }
}

int FunctionalGroup::get_total_functional_group_count(string fg_name){
    if (functional_groups->find(fg_name) == functional_groups->end() ) {
        return 0;
    } else {
        vector<FunctionalGroup*> fg = functional_groups->at(fg_name);
        int count = 0;
        for (auto fg_item : fg) {
            count += fg_item->count;
        }
        return count;
    }
}


ElementTable* FunctionalGroup::get_functional_group_elements(){
    ElementTable* _elements = create_empty_table();
    
    for (auto kv : *functional_groups){
        for (auto func_group : kv.second){
            ElementTable* fg_elements = func_group->get_elements();
            for (auto el : *fg_elements){
                _elements->at(el.first) += el.second * func_group->count;
            }
            delete fg_elements;
        }
    }
    
    return _elements;
}


void FunctionalGroup::compute_elements(){
    for (auto &kv : *functional_groups){
        for (auto func_group : kv.second){
            func_group->compute_elements();
        }
    }
}



        
string FunctionalGroup::to_string(LipidLevel level){
    string fg_string = "";
    if (is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE)){
        if ('0' <= name[0] && name[0] <= '9'){
            fg_string = (position > -1) ? (std::to_string(position) + ring_stereo + "(" + name + ")") : name;
        }
        else {
            fg_string = (position > -1) ? (std::to_string(position) + ring_stereo + name) : name;
        }
    }
    else{
        fg_string = (count > 1) ? ("(" + name + ")" + std::to_string(count)) : name;
    }
    if (stereochemistry.length() > 0 && is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE)){
        fg_string += "[" + stereochemistry + "]";
    }
            
    return fg_string;
}


int FunctionalGroup::get_double_bonds(){
    int db = count * double_bonds->get_num();
    for (auto &kv : *functional_groups){
        for (auto func_group : kv.second){
            db += func_group->get_double_bonds();
        }
    }
            
    return db;
}


void FunctionalGroup::add_position(int pos){
    position += position >= pos;
    
    for (auto &kv : *functional_groups){
        for (auto &fg : kv.second){
            fg->add_position(pos);
        }
    }
}


void FunctionalGroup::add(FunctionalGroup* fg){
    for (auto &kv : *(fg->elements)){
        elements->at(kv.first) += kv.second * fg->count;
    }
}



KnownFunctionalGroups::~KnownFunctionalGroups(){
    for (auto &kv : known_functional_groups) delete kv.second;
}

KnownFunctionalGroups KnownFunctionalGroups::k;

FunctionalGroup* KnownFunctionalGroups::get_functional_group(string fg_name){
    //static KnownFunctionalGroups k;
    if(contains_val(k.known_functional_groups, fg_name)){
        return k.known_functional_groups.at(fg_name)->copy();
    }
    return 0;
    throw RuntimeException("Name '" + fg_name + "' not registered in functional group list");
}



HeadgroupDecorator::HeadgroupDecorator(string _name, int _position, int _count, ElementTable* _elements, bool _suffix, LipidLevel _level) : FunctionalGroup(_name, _position, _count, 0, false, "", _elements){
    suffix = _suffix;
    lowest_visible_level = _level;
}

HeadgroupDecorator* HeadgroupDecorator::copy(){
    ElementTable* e = create_empty_table();
    for (auto &kv : *elements){
        e->at(kv.first) = kv.second;
    }
    return new HeadgroupDecorator(name, position, count, e, suffix, lowest_visible_level);
}



string HeadgroupDecorator::to_string(LipidLevel level){
    if (!suffix) return name + (count > 1 ? std::to_string(count) : "");

    string decorator_string = "";
    if (lowest_visible_level == NO_LEVEL || lowest_visible_level <= level){
        if (contains_val_p(functional_groups, "decorator_alkyl")){
            if (functional_groups->at("decorator_alkyl").size() > 0){
                decorator_string = (level > SPECIES) ? functional_groups->at("decorator_alkyl").at(0)->to_string(level) : "Alk";
            }
            else { 
                decorator_string = "Alk";
            }
        }
        else if (contains_val_p(functional_groups, "decorator_acyl")){
            if (functional_groups->at("decorator_acyl").size() > 0){
                decorator_string = (level > SPECIES) ? ("FA " + functional_groups->at("decorator_acyl").at(0)->to_string(level)) : "FA";
            }
            else {
                decorator_string = "FA";
            }
        }
        else {
            decorator_string = name;
        }
        decorator_string = "(" + decorator_string + ")";
    }
        
    return decorator_string;
}
