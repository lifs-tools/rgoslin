/*
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

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


#include <string>
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/LipidExceptions.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace goslin;

int main(int argc, char** argv){

    
    FattyAcid instanceZero = FattyAcid("FA1", 2, 0, 0, UNDEFINED_FA, false, 0, NULL);
    assert(0 == instanceZero.num_double_bonds);
    
    
    FattyAcid instanceOne = FattyAcid("FA1", 2, 1, 0, UNDEFINED_FA, false, 0, NULL);
    assert(1 == instanceOne.num_double_bonds);
         
    try {
        FattyAcid instanceZeroFail = FattyAcid("FA1", 2, -1, 0, UNDEFINED_FA, false, 0, NULL);
        assert(false);
    }
    catch (LipidException &e){
        assert(true);
    }
       
    FattyAcid instance = FattyAcid("FAX", 2, 0, 0, UNDEFINED_FA, false, 0, NULL);
    assert("FAX" == instance.name);

    
    instance = FattyAcid("FAX", 2, 0, 0, UNDEFINED_FA, false, 1, NULL);
    assert(1 == instance.position);
    
    try {
        instanceZero = FattyAcid("FA1", 2, 0, 0, UNDEFINED_FA, false, -2, NULL);
        assert(false);
    }
    catch (LipidException &e){
        assert(true);
    }


    instance = FattyAcid("FAX", 2, 0, 0, UNDEFINED_FA, false, 1, NULL);
    assert(2 == instance.num_carbon);


    try {
        instance = FattyAcid("FAX", 1, 0, 0, UNDEFINED_FA, false, 1, NULL);
        assert(false);
    }
    catch (LipidException &e){
        assert(true);
    }
    

    instance = FattyAcid("FAX", 2, 0, 1, UNDEFINED_FA, false, 1, NULL);
    assert(1 == instance.num_hydroxyl);


    try {
        instance = FattyAcid("FAX", 2, 0, -1, UNDEFINED_FA, false, 1, NULL);
        assert(false);
    }
    catch (LipidException &e){
        assert(true);
    }
    
    cout << "All tests passed without any problem" << endl;
    return 0;
}
