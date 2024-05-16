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


#include "cppgoslin/domain/Headgroup.h"

const map<string, vector<string> > Headgroup::glyco_table{
    {"ga1", {"Gal", "GalNAc", "Gal", "Glc"}},
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
    {"gt3", {"NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
    {"gd0a", {"HexNAc", "Hex", "NeuAc", "HexNAc", "Hex", "NeuAc", "Hex"}},
    {"gd1b", {"Gal", "GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
    {"gq1b", {"NeuAc", "NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
    {"gt1b", {"NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
    {"gd1a-ac", {"Hex", "Hex", "Hex", "HexNAc", "NeuAc", "NeuAc", "NeuAc"}},
    {"gq1-ac", {"NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
    {"gt1b-ac", {"NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}},
    {"gt3-ac", {"NeuAc", "NeuAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"}}
};
               
    
Headgroup::Headgroup(string _headgroup, vector<HeadgroupDecorator*>* _decorators, bool _use_headgroup){
    decorators = new vector<HeadgroupDecorator*>();
    
    string hg = to_lower(_headgroup);
    if (contains_val(glyco_table, hg) && !_use_headgroup){
    
        for (auto carbohydrate : glyco_table.at(hg)){
            FunctionalGroup* functional_group = 0;
            try {
                functional_group = KnownFunctionalGroups::get_functional_group(carbohydrate);
            }
            catch (const std::exception& e){
                throw LipidParsingException("Carbohydrate '" + carbohydrate + "' unknown");
            }
            
            functional_group->elements->at(ELEMENT_O) -= 1;
            decorators->push_back((HeadgroupDecorator*)functional_group);
        }
        _headgroup = "Cer";
    }
    
    
    headgroup = _headgroup;
    lipid_category = get_category(_headgroup);
    lipid_class = get_class(headgroup);
    use_headgroup = _use_headgroup;
    if (_decorators != 0){
        for (auto decorator : *_decorators) decorators->push_back(decorator);
    }
    sp_exception = lipid_category == SP && contains_val(LipidClasses::get_instance().lipid_classes.at(lipid_class).special_cases, "SP_Exception") && decorators->size() == 0;
    
}
    
Headgroup::Headgroup(Headgroup *h){
    headgroup = h->headgroup;
    lipid_category = h->lipid_category;
    lipid_class = h->lipid_class;
    use_headgroup = h->use_headgroup;
    decorators = new vector<HeadgroupDecorator*>();
    for (auto hgd : *(h->decorators)) decorators->push_back(hgd->copy());
    sp_exception = h->sp_exception;
}


Headgroup::~Headgroup(){
    for (auto hgd : *decorators) delete hgd;
    delete decorators;
}



void Headgroup::init(){
    if (!StringCategory.size()){
        for (const auto& kvp : LipidClasses::get_instance().lipid_classes){
            LipidCategory category = kvp.second.lipid_category;
            for (auto hg : kvp.second.synonyms){
                StringCategory.insert(pair<string, LipidCategory>(hg, category));
            }
        }
        
        
        for (auto kvp : LipidClasses::get_instance().lipid_classes){
            LipidClass l_class = kvp.first;
            for (auto hg : kvp.second.synonyms){
                StringClass.insert({hg, l_class});
            }
        }
        
        for (auto kvp : LipidClasses::get_instance().lipid_classes){
            ClassString.insert({kvp.first, kvp.second.synonyms.at(0)});
        }
    }
}

        

LipidCategory Headgroup::get_category(string _headgroup){
    if (!StringCategory.size()){
        for (const auto& kvp : LipidClasses::get_instance().lipid_classes){
            LipidCategory category = kvp.second.lipid_category;
            for (auto hg : kvp.second.synonyms){
                StringCategory.insert(pair<string, LipidCategory>(hg, category));
            }
        }
    }

    auto cat = StringCategory.find(_headgroup);
    return (cat != StringCategory.end()) ? StringCategory.at(_headgroup) : UNDEFINED;
}



LipidClass Headgroup::get_class(string _headgroup){
    if (!StringClass.size()){
        for (auto kvp : LipidClasses::get_instance().lipid_classes){
            LipidClass l_class = kvp.first;
            for (auto hg : kvp.second.synonyms){
                StringClass.insert({hg, l_class});
            }
        }
    }
    
    auto cl = StringClass.find(_headgroup);
    return (cl != StringClass.end()) ? cl->second : UNDEFINED_CLASS;
}


string Headgroup::get_class_string(LipidClass _lipid_class){
    if (!ClassString.size()){
        for (auto kvp : LipidClasses::get_instance().lipid_classes){
            ClassString.insert({kvp.first, kvp.second.synonyms.at(0)});
        }
    }
    auto cl = ClassString.find(_lipid_class);
    return (cl != ClassString.end()) ? ClassString.at(_lipid_class) : "UNDEFINED";
}


string Headgroup::get_class_name(){
    if (uncontains_val(LipidClasses::get_instance().lipid_classes, lipid_class)) return "UNDEFINED";
    return LipidClasses::get_instance().lipid_classes.at(lipid_class).class_name;
}


string Headgroup::get_category_string(LipidCategory _lipid_category){
    return CategoryString.at(_lipid_category);
}


bool Headgroup::decorator_sorting (HeadgroupDecorator* hi, HeadgroupDecorator* hj){
    return hi->name < hj->name;
}
        
        
string Headgroup::get_lipid_string(LipidLevel level){
    if (level == CATEGORY){
        return get_category_string(lipid_category);
    }
    
    string hgs = ((use_headgroup) ? headgroup : get_class_string(lipid_class));
    
    
    stringstream headgoup_string;
         
    // adding prefixes to the headgroup
    if (!is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE | STRUCTURE_DEFINED)){
        vector<HeadgroupDecorator*> decorators_sorted;
        for (auto hgd : *decorators){
            if (hgd->suffix) continue;
            
            HeadgroupDecorator* hgd_copy = hgd->copy();
            hgd_copy->name = goslin::replace_all(hgd_copy->name, "Gal", "Hex");
            hgd_copy->name = goslin::replace_all(hgd_copy->name, "Glc", "Hex");
            hgd_copy->name = goslin::replace_all(hgd_copy->name, "S(3')", "S");
            
            decorators_sorted.push_back(hgd_copy);
        }
        
        if (decorators_sorted.size() > 0){
            sort (decorators_sorted.begin(), decorators_sorted.end(), decorator_sorting);
            
            for (int i = decorators_sorted.size() - 1; i > 0; --i){
                HeadgroupDecorator* hge = decorators_sorted[i];
                HeadgroupDecorator* hge_before = decorators_sorted[i - 1];
                if (hge->name == hge_before->name){
                    hge_before->count += hge->count;
                    delete hge;
                    decorators_sorted.erase(decorators_sorted.begin() + i);
                }
            }
            for (auto hge : decorators_sorted){
                headgoup_string << hge->to_string(level);
                delete hge;
            }
        }
    }
    else {
        for (auto hgd : *decorators){
            if (!hgd->suffix) headgoup_string << hgd->to_string(level) << "-";
        }
    }
    
    // adding headgroup
    headgoup_string << hgs;
    // ading suffixes to the headgroup
    for (auto hgd : *decorators){
        if (hgd != 0 && hgd->suffix){
            headgoup_string << hgd->to_string(level);
        }
    }
    if (is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE) && lipid_category == SP && !sp_exception){
        headgoup_string << "(1)";
    }
    
    return headgoup_string.str();
}
        
        
ElementTable* Headgroup::get_elements(){
    ClassMap &lipid_classes = LipidClasses::get_instance().lipid_classes;
    
    if (use_headgroup || uncontains_val(lipid_classes, lipid_class)){
        throw RuntimeException("Element table cannot be computed for lipid '" + headgroup + "'");
    }
    
    ElementTable *elements = create_empty_table();
    
    for (auto &kv : lipid_classes.at(lipid_class).elements){
        elements->at(kv.first) += kv.second;
    }
    
    
    for (auto hgd : *decorators){
        ElementTable *hgd_elements = hgd->get_elements();
        for (auto &kv : *hgd_elements){
            elements->at(kv.first) += kv.second * hgd->count;
        }
        delete hgd_elements;
    }
    
    return elements;
}
