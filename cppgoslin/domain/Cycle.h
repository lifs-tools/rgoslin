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


#ifndef CYCLE_H
#define CYCLE_H

#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/Element.h"
#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/DoubleBonds.h"
#include "cppgoslin/domain/StringFunctions.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cctype>
#include <string>

using namespace std;
using namespace goslin;

class Cycle : public FunctionalGroup {
public:
    int cycle;
    int start;
    int end;
    vector<Element>* bridge_chain;
    
    Cycle(int _cycle, int _start = -1, int _end = -1, DoubleBonds* _double_bonds = 0, map<string, vector<FunctionalGroup*> >* _functional_groups = 0, vector<Element>* _bridgeChain = 0);
    ~Cycle();
    Cycle* copy();
    int get_double_bonds();
    void rearrange_functional_groups(FunctionalGroup *parent, int shift);
    void shift_positions(int shift);
    void compute_elements();
    void add_position(int pos);
    string to_string(LipidLevel level);
};

#endif /* CYCLE_H */
