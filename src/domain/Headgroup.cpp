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
    
Headgroup::Headgroup(string _headgroup, vector<HeadgroupDecorator*>* _decorators, bool _use_headgroup){
    headgroup = _headgroup;
    lipid_category = get_category(_headgroup);
    lipid_class = get_class(headgroup);
    use_headgroup = _use_headgroup;
    decorators = new vector<HeadgroupDecorator*>();
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
    
    if (level == CLASS){
        return hgs;
    }
    
    stringstream headgoup_string;
         
    // adding prefixes to the headgroup
    if (!is_level(level, COMPLETE_STRUCTURE | FULL_STRUCTURE | STRUCTURE_DEFINED)){
        vector<HeadgroupDecorator*> decorators_tmp;
        for (auto hgd : *decorators){
            if (!hgd->suffix) decorators_tmp.push_back(hgd->copy());
        }
        sort (decorators_tmp.begin(), decorators_tmp.end(), decorator_sorting);
        
        for (int i = decorators_tmp.size() - 1; i > 0; --i){
            HeadgroupDecorator* hge = decorators_tmp[i];
            HeadgroupDecorator* hge_before = decorators_tmp[i - 1];
            if (hge->name == hge_before->name){
                hge_before->count += hge->count;
                delete hge;
                decorators_tmp.erase(decorators_tmp.begin() + i);
            }
        }
        for (auto hge : decorators_tmp){
            headgoup_string << hge->to_string(level);
            delete hge;
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
