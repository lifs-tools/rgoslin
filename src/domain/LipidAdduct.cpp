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


#include "cppgoslin/domain/LipidAdduct.h"


LipidAdduct::LipidAdduct(){
    lipid = NULL;
    adduct = NULL;
    sum_formula = "";
}

LipidAdduct::~LipidAdduct(){
    if (lipid) delete lipid;
    if (adduct) delete adduct;
}
    
    
string LipidAdduct::get_lipid_string(LipidLevel level){
    stringstream s;
    if (lipid) s << lipid->get_lipid_string(level);
    else return "";
    
    
    switch (level){
        case CLASS:
        case CATEGORY:
            break;
            
        default:
            if (adduct) s << adduct->get_lipid_string();
            break;
    }
    
    return s.str();
}


string LipidAdduct::get_class_name(){
    return (lipid) ? lipid->headgroup->get_class_name() : "";
}



LipidLevel LipidAdduct::get_lipid_level(){
    return lipid->get_lipid_level();
}


    
string LipidAdduct::get_extended_class(){
    return lipid ? lipid->get_extended_class() : "";
}



void LipidAdduct::sort_fatty_acyl_chains(){
    lipid->sort_fatty_acyl_chains();
}



double LipidAdduct::get_mass(){
    ElementTable* elements = get_elements();
    int charge = (adduct != NULL) ? adduct->get_charge() : 0;
    double mass = goslin::get_mass(elements);
    if (charge != 0) mass = (mass - charge * ELECTRON_REST_MASS) / fabs(charge);
    delete elements;
    
    return mass;
}


ElementTable* LipidAdduct::get_elements(){
    ElementTable* elements = create_empty_table();
    
    if (lipid != NULL){
        ElementTable* lipid_elements = lipid->get_elements();
        for (auto e : *lipid_elements) elements->at(e.first) += e.second;
        
        delete lipid_elements;
    }
            
    if (adduct != NULL){
        ElementTable* adduct_elements = adduct->get_elements();
        for (auto e : *adduct_elements) elements->at(e.first) += e.second;
            
        delete adduct_elements;
    }
    return elements;
}
    
    
string LipidAdduct::get_lipid_fragment_string(LipidLevel level){
    stringstream s;
    
    if (lipid) s << lipid->get_lipid_string(level);
    else return "";
    
    
    
    switch (level){
        case CLASS:
        case CATEGORY:
            break;
    
        default:
            if (adduct) s << adduct->get_lipid_string();
            break;
    }
    
    return s.str();
}


string LipidAdduct::get_sum_formula(){
    ElementTable* elements = get_elements();
    
    string formula = compute_sum_formula(elements);
    delete elements;
            
    return formula;
}
