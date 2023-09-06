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


#ifndef PARSER_CLASSES_H
#define PARSER_CLASSES_H


#include "cppgoslin/domain/StringFunctions.h"
#include <string>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>


using namespace std;

template<class T>
class BaseParserEventHandler;


class GrammarString : public string {
public:
    GrammarString(string s) : string(s){}
        
};



    
// DP stands for dynamic programming
class DPNode {
public:    
    uint64_t rule_index_1;
    uint64_t rule_index_2;
    DPNode *left = NULL;
    DPNode *right = NULL;
    
    DPNode(uint64_t _rule1, uint64_t _rule2, DPNode *_left, DPNode *_right);
};
        

class TreeNode {
public:
    uint64_t rule_index;
    TreeNode *left;
    TreeNode *right;
    char terminal;
    bool fire_event;
    static const char EOF_SIGN = '\0';
    
    TreeNode(uint64_t _rule, bool _fire_event);
    ~TreeNode();
    string get_text();
    int get_int();
};
   

        

// this class is dedicated to have an efficient sorted set class storing
// values within 0..n-1 and fast sequencial iterator
class Bitfield {
    class iter;
    
public:
    uint64_t *field;
    uint64_t field_len;
    uint64_t num_size;
    uint64_t length;
    
    iter begin();
    iter end();
    
    uint64_t size() const;
    
    
    Bitfield(uint64_t length, bool filled_with_ones = false);
    ~Bitfield();
    void insert(uint64_t pos);
    void remove(uint64_t pos);
    bool find(uint64_t pos);
    void init();
    int next(int pos = -1);
    void print_bitfield(uint64_t l);
    uint64_t size();
    
private:
    class iter : public std::iterator<std::output_iterator_tag, int>{
        public:
            explicit iter(Bitfield& _bitfield, uint64_t index = 0);
            int operator*();
            iter & operator++();
            iter & operator++(int);
            bool operator!=(const iter &) const;
            uint64_t num_index;
            int last_position;
            Bitfield &bitfield;
            bool get_next;
    };
};


#endif /* PARSER_CLASSES_H */


