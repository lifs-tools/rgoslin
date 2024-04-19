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


#include "cppgoslin/domain/GenericDatastructures.h"


GenericList::~GenericList(){
    for (auto& x : list){
        del(x);
    }
}



void GenericList::remove_all(){
    while(!list.empty()){
        del(list.back());
        list.pop_back();
    }
}


void GenericList::del(pair<int, void*> &x){
    switch (x.first){
        case TYPE_BOOL:
            delete (bool*)x.second;
            break;
            
        case TYPE_INT:
            delete (int*)x.second;
            break;
            
        case TYPE_LONG:
            delete (long*)x.second;
            break;
            
        case TYPE_FLOAT:
            delete (float*)x.second;
            break;
            
        case TYPE_DOUBLE:
            delete (double*)x.second;
            break;
            
        case TYPE_STRING:
            delete (string*)x.second;
            break;
            
        case TYPE_LIST:
            delete (GenericList*)x.second;
            break;
            
        case TYPE_DICTIONARY:
            delete (GenericList*)x.second;
            break;
            
        default:
            break;
    }
}



void GenericList::add_int(int i){
    void* vi = new int(i);
    list.push_back({TYPE_INT, vi});
}

void GenericList::set_int(int i, int ii){
    void* vi = new int(ii);
    del(list.at(i));
    list.at(i) = {TYPE_INT, vi};
}

int GenericList::get_int(int key){
    key += list.size() * (key < 0);
    return *((int*)(list.at(key).second));
}

void GenericList::add_long(long l){
    void* vl = new long(l);
    list.push_back({TYPE_LONG, vl});
}

void GenericList::set_long(int i, long l){
    void* vl = new long(l);
    del(list.at(i));
    list.at(i) = {TYPE_LONG, vl};
}

long GenericList::get_long(int key){
    key += list.size() * (key < 0);
    return *((long*)(list.at(key).second));
}

void GenericList::add_float(float f){
    void* vf = new float(f);
    list.push_back({TYPE_FLOAT, vf});
}

void GenericList::set_float(int i, float f){
    void* vf = new float(f);
    del(list.at(i));
    list.at(i) = {TYPE_FLOAT, vf};
}

float GenericList::get_float(int key){
    key += list.size() * (key < 0);
    return *((float*)(list.at(key).second));
}

void GenericList::add_double(double d){
    void* vd = new double(d);
    list.push_back({TYPE_DOUBLE, vd});
}

void GenericList::set_double(int i, double d){
    void* vd = new double(d);
    del(list.at(i));
    list.at(i) = {TYPE_DOUBLE, vd};
}

double GenericList::get_double(int key){
    key += list.size() * (key < 0);
    return *((double*)(list.at(key).second));
}

void GenericList::add_string(string s){
    void* vs = new string(s);
    list.push_back({TYPE_STRING, vs});
}

void GenericList::set_string(int i, string s){
    void* vs = new string(s);
    del(list.at(i));
    list.at(i) = {TYPE_STRING, vs};
}

string GenericList::get_string(int key){
    key += list.size() * (key < 0);
    return *((string*)(list.at(key).second));
}

void GenericList::add_list(GenericList* v){
    list.push_back({TYPE_LIST, v});
}

void GenericList::set_list(int i, GenericList* v){
    del(list.at(i));
    list.at(i) = {TYPE_LIST, v};
}

GenericList* GenericList::get_list(int key){
    key += list.size() * (key < 0);
    return (GenericList*)(list.at(key).second);
}


void GenericList::add_dictionary(GenericDictionary* dict){
    list.push_back({TYPE_DICTIONARY, (void*)dict});
}


void GenericList::set_dictionary(int i, GenericDictionary* dict){
    del(list.at(i));
    list.at(i) = {TYPE_DICTIONARY, (void*)dict};
}


GenericDictionary* GenericList::get_dictionary(int key){
    key += list.size() * (key < 0);
    return (GenericDictionary*)(list.at(key).second);
}






GenericDictionary::~GenericDictionary(){
    for (auto& x : dictionary){
        del(x.second);
    }
}

void GenericDictionary::remove(string key){
    del(dictionary.at(key));
    dictionary.erase(key);
}



void GenericDictionary::remove_all(){
    vector<string> keys;
    for (auto &kv : dictionary) keys.push_back(kv.first);
    for (auto key : keys){
        del(dictionary.at(key));
        dictionary.erase(key);
    }
}


void GenericDictionary::del(pair<int, void*> &x){
    switch (x.first){
        case TYPE_BOOL:
            delete (bool*)x.second;
            break;
            
        case TYPE_INT:
            delete (int*)x.second;
            break;
            
        case TYPE_LONG:
            delete (long*)x.second;
            break;
            
        case TYPE_FLOAT:
            delete (float*)x.second;
            break;
            
        case TYPE_DOUBLE:
            delete (double*)x.second;
            break;
            
        case TYPE_STRING:
            delete (string*)x.second;
            break;
            
        case TYPE_LIST:
            delete (GenericList*)x.second;
            break;
            
        case TYPE_DICTIONARY:
            delete (GenericDictionary*)x.second;
            break;
            
        default:
            break;
    }
}

void GenericDictionary::set_int(string key, int i){
    void* vi = new int(i);
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, vi};
    }
    else {
        dictionary.insert({key, {TYPE_INT, vi}});
    }
}

int GenericDictionary::get_int(string key){
    return *((int*)(dictionary.at(key).second));
}

void GenericDictionary::set_long(string key, long l){
    void* vl = new long(l);
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, vl};
    }
    else {
        dictionary.insert({key, {TYPE_LONG, vl}});
    }
}

long GenericDictionary::get_long(string key){
    return *((long*)(dictionary.at(key).second));
}

void GenericDictionary::set_float(string key, float f){
    void* vf = new float(f);
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, vf};
    }
    else {
        dictionary.insert({key, {TYPE_FLOAT, vf}});
    }
}

float GenericDictionary::get_float(string key){
    return *((float*)(dictionary.at(key).second));
}

void GenericDictionary::set_double(string key, double d){
    void* vd = new double(d);
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, vd};
    }
    else {
        dictionary.insert({key, {TYPE_DOUBLE, vd}});
    }
}

double GenericDictionary::get_double(string key){
    return *((double*)(dictionary.at(key).second));
}

void GenericDictionary::set_string(string key, string s){
    void* vs = new string(s);
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, vs};
    }
    else {
        dictionary.insert({key, {TYPE_STRING, vs}});
    }
}

string GenericDictionary::get_string(string key){
    return *((string*)(dictionary.at(key).second));
}

void GenericDictionary::set_list(string key, GenericList* v){
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, v};
    }
    else {
        dictionary.insert({key, {TYPE_LIST, v}});
    }
}

GenericList* GenericDictionary::get_list(string key){
    return (GenericList*)(dictionary.at(key).second);
}


void GenericDictionary::set_dictionary(string key, GenericDictionary* dict){
    if (contains_val(dictionary, key)){
        del(dictionary.at(key));
        dictionary.at(key) = {TYPE_INT, dict};
    }
    else {
        dictionary.insert({key, {TYPE_DICTIONARY, (void*)dict}});
    }
}


GenericDictionary* GenericDictionary::get_dictionary(string key){
    return (GenericDictionary*)(dictionary.at(key).second);
}

bool GenericDictionary::contains_key(string key){
    return contains_val(dictionary, key);
}

