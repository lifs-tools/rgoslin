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


#ifndef SUM_FORMULA_PARSER_EVENT_HANDLER_H
#define SUM_FORMULA_PARSER_EVENT_HANDLER_H

#include "cppgoslin/domain/Element.h"
#include "cppgoslin/parser/BaseParserEventHandler.h"
#include <string>
#include <set>
#include <map>
#include <vector>

using namespace std;
using namespace goslin;

class SumFormulaParserEventHandler : public BaseParserEventHandler<ElementTable*> {
public:
    Element element;
    int count;
    
    SumFormulaParserEventHandler();
    void reset_parser(TreeNode *node);
    void element_group_post_event(TreeNode *node);
    void element_pre_event(TreeNode *node);
    void single_element_group_pre_event(TreeNode *node);
    void count_pre_event(TreeNode *node);
};


#endif /* SUM_FORMULA_PARSER_EVENT_HANDLER_H */
        
