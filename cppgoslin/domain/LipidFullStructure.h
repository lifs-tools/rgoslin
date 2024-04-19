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


#ifndef LIPID_FULL_STRUCTURE_H
#define LIPID_FULL_STRUCTURE_H

#include <string>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidStructureDefined.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/Headgroup.h"
#include <sstream>
#include <vector>

using namespace std;
using namespace goslin;

class LipidFullStructure : public LipidStructureDefined {
public:
    LipidFullStructure(Headgroup* _headgroup, vector<FattyAcid*>* _fa = NULL);
    string get_lipid_string(LipidLevel level = NO_LEVEL);
    LipidLevel get_lipid_level();
            
};

#endif /* LIPID_FULL_STRUCTURE_H */
