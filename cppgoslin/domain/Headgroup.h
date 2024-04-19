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


#ifndef HEADGROUP_H
#define HEADGROUP_H

#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/LipidSpeciesInfo.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/domain/Element.h"
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

class Headgroup {
public:
    string headgroup;
    LipidCategory lipid_category;
    LipidClass lipid_class;
    bool use_headgroup;
    vector<HeadgroupDecorator*>* decorators;
    bool sp_exception;
    const set<string> exception_headgroups {"Cer", "SPB"};
    static const map<string, vector<string> > glyco_table;
    
    Headgroup(string _headgroup, vector<HeadgroupDecorator*>* _decorators = 0, bool _use_headgroup = false);
    Headgroup(Headgroup *h);
    ~Headgroup();
    static void init();
    string get_lipid_string(LipidLevel level = NO_LEVEL);
    ElementTable* get_elements();
    static LipidCategory get_category(string _headgroup);
    static LipidClass get_class(string _head_group);
    static string get_class_string(LipidClass _lipid_class);
    static string get_category_string(LipidCategory _lipid_category);
    static bool decorator_sorting (HeadgroupDecorator* hi, HeadgroupDecorator* hj);
    string get_class_name();
    
};

#endif /* HEADGROUP_H */
