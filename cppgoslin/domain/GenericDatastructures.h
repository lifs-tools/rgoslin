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


#ifndef GENERICDATASTRUCTURES_H
#define GENERICDATASTRUCTURES_H

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <tuple>
#include "cppgoslin/domain/StringFunctions.h"

enum Type {TYPE_INT, TYPE_LONG, TYPE_FLOAT, TYPE_DOUBLE, TYPE_BOOL, TYPE_STRING, TYPE_LIST, TYPE_DICTIONARY};

using namespace std;

class GenericDictionary;

class GenericList {
public:
    ~GenericList();
    vector<pair<int, void*>> list;
    
    void add_int(int i);
    void set_int(int i, int ii);
    int get_int(int i);
    
    void add_long(long l);
    void set_long(int i, long l);
    long get_long(int i);
    
    void add_float(float f);
    void set_float(int i, float f);
    float get_float(int i);
    
    void add_double(double d);
    void set_double(int i, double d);
    double get_double(int i);
    
    void add_string(string s);
    void set_string(int i, string s);
    string get_string(int i);
    
    void add_list(GenericList* v);
    void set_list(int i, GenericList* v);
    GenericList* get_list(int i);
    
    void add_dictionary(GenericDictionary* d);
    void set_dictionary(int i, GenericDictionary* d);
    GenericDictionary* get_dictionary(int i);
    
    void del(pair<int, void*> &x);
    void remove_all();
};




class GenericDictionary {
public:
    ~GenericDictionary();
    map<string, pair<int, void*>> dictionary;
    
    void set_int(string key, int i);
    int get_int(string key);
    
    void set_long(string key, long l);
    long get_long(string key);
    
    void set_float(string key, float f);
    float get_float(string key);
    
    void set_double(string key, double d);
    double get_double(string key);
    
    void set_string(string key, string s);
    string get_string(string key);
    
    void set_list(string key, GenericList* v);
    GenericList* get_list(string key);
    
    void set_dictionary(string key, GenericDictionary* dict);
    GenericDictionary* get_dictionary(string key);
    
    void remove(string key);
    void remove_all();
    void del(pair<int, void*> &x);
    bool contains_key(string key);
};

#endif /* GENERICDATASTRUCTURES_H */
