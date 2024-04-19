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


#ifndef ADDUCT_H
#define ADDUCT_H

#include <string>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/Element.h"
#include "cppgoslin/domain/StringFunctions.h"
#include <sstream>
#include <map>


using namespace goslin;
using namespace std;





class Adduct {
public:
    string sum_formula;
    string adduct_string;
    int charge;
    int charge_sign;
    ElementTable heavy_elements;
    
    static const map<string, int> adduct_charges;
    
    Adduct(string _sum_formula, string _adduct_string, int _charge = 0, int _sign = 1);
    Adduct(Adduct *a);
    void set_charge_sign(int sign);
    string get_lipid_string();
    ElementTable* get_elements();
    int get_charge();
    string get_heavy_isotope_string();
};



class KnownAdducts {
    public:
        static KnownAdducts& get_instance()
        {
            static KnownAdducts instance;
            return instance;
        }
    private:
        KnownAdducts();
        
    public:
        map<string, ElementTable> known_adducts;
        KnownAdducts(KnownAdducts const&) = delete;
        void operator=(KnownAdducts const&) = delete;
};


#endif /* ADDUCT_H */
