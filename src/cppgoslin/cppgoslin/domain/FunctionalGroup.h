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


#ifndef FUNCTIONALGROUP_H
#define FUNCTIONALGROUP_H

#include "cppgoslin/domain/Element.h"
#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/DoubleBonds.h"
#include "cppgoslin/domain/StringFunctions.h"
#include <map>
#include <vector>
#include <string>


using namespace std;
using namespace goslin;

class FunctionalGroup {
public:
    string name;
    int position;
    int count;
    int num_atoms;
    string stereochemistry;
    string ring_stereo;
    DoubleBonds* double_bonds;
    bool is_atomic;
    ElementTable* elements;
    map<string, vector<FunctionalGroup*>>* functional_groups;
    static bool position_sort_function(FunctionalGroup* f1, FunctionalGroup *f2);
    static bool lower_name_sort_function(string s1, string s2);
    
    FunctionalGroup(string _name, int _position = -1, int _count = 1, DoubleBonds* _double_bonds = 0, bool _is_atomic = false, string _stereochemistry = "", ElementTable* _elements = 0, map<string, vector<FunctionalGroup*>>* _functional_groups = 0);
    virtual ~FunctionalGroup();
    virtual FunctionalGroup* copy();
    virtual ElementTable* get_elements();
    virtual void shift_positions(int shift);
    virtual int get_total_functional_group_count(string fg_name);
    virtual ElementTable* get_functional_group_elements();
    virtual void compute_elements();
    virtual string to_string(LipidLevel level);
    virtual int get_double_bonds();
    virtual void add_position(int pos);
    void add(FunctionalGroup* fg);
};


class KnownFunctionalGroups {
public:
    KnownFunctionalGroups();
    ~KnownFunctionalGroups();
    map<string, FunctionalGroup*> known_functional_groups;
    static FunctionalGroup* get_functional_group(string fg_name);
    static KnownFunctionalGroups k;
};


class HeadgroupDecorator : public FunctionalGroup {
public:
    bool suffix;
    LipidLevel lowest_visible_level;
    
    HeadgroupDecorator(string _name, int _position = -1, int _count = 1, ElementTable* _elements = 0, bool _suffix = false, LipidLevel _level = NO_LEVEL);
    string to_string(LipidLevel level);
    HeadgroupDecorator* copy();
};

    

#endif /* FUNCTIONALGROUP_H */
        
        
