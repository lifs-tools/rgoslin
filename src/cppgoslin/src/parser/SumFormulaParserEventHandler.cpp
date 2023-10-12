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


#include "cppgoslin/parser/SumFormulaParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<ElementTable*>::registered_events->insert({x, bind(&SumFormulaParserEventHandler::y, this, placeholders::_1)})
    

SumFormulaParserEventHandler::SumFormulaParserEventHandler() : BaseParserEventHandler<ElementTable*>() {
    content = create_empty_table();
    element = ELEMENT_H;
    count = 0;
    
    reg("molecule_pre_event", reset_parser);
    reg("element_group_post_event", element_group_post_event);
    reg("element_pre_event", element_pre_event);
    reg("single_element_pre_event", single_element_group_pre_event);
    reg("count_pre_event", count_pre_event);
    
}


void SumFormulaParserEventHandler::reset_parser(TreeNode *node){
    content = create_empty_table();
}


void SumFormulaParserEventHandler::element_group_post_event(TreeNode *node){
    content->at(element) += count;
}


void SumFormulaParserEventHandler::element_pre_event(TreeNode *node){
    string parsed_element = node->get_text();
    
    if (element_positions.find(parsed_element) != element_positions.end()){
        element = element_positions.at(parsed_element);
    }
            
    else {
        throw LipidException("Error: element '" + parsed_element + "' is unknown");
    }
}


void SumFormulaParserEventHandler::single_element_group_pre_event(TreeNode *node){
    string parsed_element = node->get_text();
    if (element_positions.find(parsed_element) != element_positions.end()){
        element = element_positions.at(parsed_element);
        content->at(element) += 1;
    }
        
    else {
        throw LipidException("Error: element '" + parsed_element + "' is unknown");
    }
}


void SumFormulaParserEventHandler::count_pre_event(TreeNode *node){
    count = node->get_int();
}


