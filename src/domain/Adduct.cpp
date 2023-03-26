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


#include "cppgoslin/domain/Adduct.h"



KnownAdducts::KnownAdducts(){
    known_adducts = {
        {"+H", {{ELEMENT_H, 1}} },
        {"+2H", {{ELEMENT_H, 2}} },
        {"+3H", {{ELEMENT_H, 3}} },
        {"+4H", {{ELEMENT_H, 4}} },
        {"-H", {{ELEMENT_H, -1}} },
        {"-2H", {{ELEMENT_H, -2}} },
        {"-3H", {{ELEMENT_H, -3}} },
        {"-4H", {{ELEMENT_H, -4}} },
        {"+H-H2O", {{ELEMENT_H, -1}, {ELEMENT_O, -1}} },
        {"+NH4", {{ELEMENT_N, 1}, {ELEMENT_H, 4}} },
        {"+Cl", {{ELEMENT_Cl, 1}} },
        {"+HCOO", {{ELEMENT_H, 1}, {ELEMENT_C, 1}, {ELEMENT_O, 2}} },
        {"+CH3COO", {{ELEMENT_H, 3}, {ELEMENT_C, 2}, {ELEMENT_O, 2}} },
    };
}


Adduct::Adduct(Adduct *a){
    if (a != 0){
        sum_formula = a->sum_formula;
        adduct_string = a->adduct_string;
        charge = a->charge;
        charge_sign = a->charge_sign;
        for (auto e : element_order) heavy_elements.insert({e, a->heavy_elements.at(e)});
    }
}

    
Adduct::Adduct(string _sum_formula, string _adduct_string, int _charge, int _sign){
    sum_formula = _sum_formula;
    adduct_string = _adduct_string;
    charge = _charge;
    set_charge_sign(_sign);
    for (auto &e : element_order) heavy_elements.insert({e, 0});
}

const map<string, int> Adduct::adduct_charges {
    {"+H", 1},  {"+2H", 2}, {"+3H", 3}, {"+4H", 4},
    {"-H", -1}, {"-2H", -2}, {"-3H", -3}, {"-4H", -4},
    {"+H-H2O", 1}, {"+NH4", 1}, {"+Cl", -1}, {"+HCOO", -1}, {"+CH3COO", -1}
};


string Adduct::get_heavy_isotope_string(){
    stringstream ss;
    for (auto e : element_order){
        if (heavy_elements[e] > 0){
            ss << heavy_elements[e] << heavy_shortcut.at(e);
        }
    }
    return ss.str();
}

void Adduct::set_charge_sign(int sign){
    if (-1 <= sign && sign <= 1){
        charge_sign = sign;
    }
        
    else {
        throw ConstraintViolationException("Sign can only be -1, 0, or 1");
    }
}
        
string Adduct::get_lipid_string(){
    if (charge == 0){
        return "[M" + get_heavy_isotope_string() + "]";
    }
    stringstream stst;
    stst << "[M" << get_heavy_isotope_string() << sum_formula << adduct_string << "]" << charge << ((charge_sign > 0) ? "+" : "-");
    
    return stst.str();
}

ElementTable* Adduct::get_elements(){
    ElementTable* elements = create_empty_table();
    
    for (auto &kv : heavy_elements){
        if (kv.second > 0){
            elements->at(heavy_to_regular.at(kv.first)) -= kv.second;
            elements->at(kv.first) += kv.second;
        }
    }
    
    if (adduct_string.length() > 0){
        if (contains_val(adduct_charges, adduct_string)){
            if (adduct_charges.at(adduct_string) != get_charge()){
                throw ConstraintViolationException("Provided charge '" + std::to_string(get_charge()) + "' in contradiction to adduct '" + adduct_string + "' charge '" + std::to_string(adduct_charges.at(adduct_string)) + "'.");
            }
            for (auto kv : KnownAdducts::get_instance().known_adducts.at(adduct_string)){
                elements->at(kv.first) += kv.second;
            }
        }
        else {
            throw ConstraintViolationException("Adduct '" + adduct_string + "' is unknown.");
        }
    }
    
    return elements;
}



int Adduct::get_charge(){
    return charge * charge_sign;
}



