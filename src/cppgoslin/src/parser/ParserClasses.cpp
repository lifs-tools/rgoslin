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


#include "cppgoslin/parser/Parser.h"



DPNode::DPNode(uint64_t _rule1, uint64_t _rule2, DPNode *_left, DPNode *_right){
    rule_index_1 = _rule1;
    rule_index_2 = _rule2;
    left = _left;
    right = _right;
}

      
      
      
    
TreeNode::TreeNode(uint64_t _rule, bool _fire_event){
    rule_index = _rule;
    left = NULL;
    right = NULL;
    terminal = 0;
    fire_event = _fire_event;
}





TreeNode::~TreeNode(){
    if (left != NULL) delete left;
    if (right != NULL) delete right;
}




string TreeNode::get_text(){
    if (!terminal){
        string left_str = left->get_text();
        string right_str = right != NULL ? right->get_text() : "";
        return (left_str != string(1, EOF_SIGN) ? left_str : "") + (right_str != string(1, EOF_SIGN) ? right_str : "");
    }
    return string(1, (char)terminal);
}

int TreeNode::get_int(){
    if (!terminal){
        string left_str = left->get_text();
        string right_str = right != NULL ? right->get_text() : "";
        return atoi(((left_str != string(1, EOF_SIGN) ? left_str : "") + (right_str != string(1, EOF_SIGN) ? right_str : "")).c_str());
    }
    return atoi(string(1, (char)terminal).c_str());
}
       
       
       
       

Bitfield::Bitfield(uint64_t _length, bool filled_with_ones){
    length = _length;
    field_len = 1 + ((length + 1) >> 6);
    field = new uint64_t[field_len];
    num_size = 0;
    if (filled_with_ones){
        num_size = (field_len << 6) - 1;
        for (uint64_t i = 0; i < field_len; ++i) field[i] = ~(0ull);
        
        while (num_size > length){
            field[num_size >> 6] &= ~((uint64_t)(1ull << (num_size & 63)));
            --num_size;
        }
    }
    else {
        num_size = 0;
        for (uint64_t i = 0; i < field_len; ++i) field[i] = 0ull;
    }
}



Bitfield::~Bitfield(){
    delete []field;
}




void Bitfield::insert(uint64_t pos){
    if (!find(pos)){
        field[pos >> 6] |= (uint64_t)(1ull << (pos & 63));
        ++num_size;
    }
}


void Bitfield::remove(uint64_t pos){
    if (find(pos)){
        field[pos >> 6] &= ~((uint64_t)(1ull << (pos & 63)));
        --num_size;
    }
}




bool Bitfield::find(uint64_t pos){
    if (pos > length) return false;
    return ((field[pos >> 6] >> (pos & 63)) & 1ull) == 1ull;
}




void Bitfield::print_bitfield(uint64_t l){
    for (int i = 63; i >= 0; --i){
        cout << ((l >> i) & 1);
    } cout << endl;
}



int Bitfield::next(int pos){
    if ((int)pos >= (int)length) throw RuntimeException("Bitfield out of range");
    
    
    uint64_t field_pos = pos >> 6;
    uint64_t field_bits = field[field_pos] & (~((1ull << (pos & 63)) - 1ull));
    
    do {
        if (field_bits){
            return (field_pos << 6) + __builtin_ctzll(field_bits & -field_bits);
        }
        if (++field_pos < field_len) field_bits = field[field_pos];
    } while (field_pos < field_len);
    
    throw RuntimeException("Bitfield out of range");
}





Bitfield::iter::iter(Bitfield & _bitfield, uint64_t index) : bitfield (_bitfield) {
    num_index = index;
    last_position = -1;
    get_next = true;
}




int Bitfield::iter::operator*() {
    if (get_next){
        last_position = bitfield.next(last_position + 1);
        get_next = false;
    }
    return last_position;
}




Bitfield::iter & Bitfield::iter::operator++() {
    num_index++;
    get_next = true;
    return *this;
}




Bitfield::iter & Bitfield::iter::operator++(int i) {
    return ++(*this);
}




bool Bitfield::iter::operator!=(const iter & rhs) const {
    return num_index != rhs.num_index;
}




uint64_t Bitfield::size() {
    return num_size;
}




Bitfield::iter Bitfield::begin() {
    return Bitfield::iter(*this, 0);
}





Bitfield::iter Bitfield::end(){
    return Bitfield::iter(*this, size());
}


