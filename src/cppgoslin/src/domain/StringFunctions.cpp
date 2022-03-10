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


#include "cppgoslin/domain/StringFunctions.h"

using namespace goslin;

string goslin::compute_sum_formula(ElementTable* elements){
    stringstream ss;
    
    for (auto e : element_order){
        if (elements->at(e) > 0) ss << element_shortcut.at(e);
        if (elements->at(e) > 1) ss << elements->at(e);
    }
    return ss.str();
}


string goslin::to_lower(string st){
    string s = string(st);
    std::transform(s.begin(), s.end(), s.begin(),[](unsigned char c){ return std::tolower(c); });
    return s;
}


string goslin::to_upper(string st){
    string s = string(st);
    std::transform(s.begin(), s.end(), s.begin(),[](unsigned char c){ return std::toupper(c); });
    return s;
}

bool goslin::endswith(const string &main_str, const string &to_match){
    if(main_str.size() >= to_match.size() &&
            main_str.compare(main_str.size() - to_match.size(), to_match.size(), to_match) == 0)
            return true;
        else
            return false;
}


string goslin::strip(string s, char c){
    if (s.length() > 0) {
        uint32_t st = 0;
        while (st < s.length() - 1 && s[st] == c) ++st;
        s = s.substr(st, s.length() - st);
    }
    
    if (s.length() > 0) {
        uint32_t en = 0;
        while (en < s.length() - 1 && s[s.length() - 1 - en] == c) ++en;
        s = s.substr(0, s.length() - en);
    }
    return s;
}



double goslin::get_mass(ElementTable *elements){
    double mass = 0;
    for (auto e : *elements) mass += element_masses.at(e.first) * e.second;
    return mass;
}



ElementTable* goslin::create_empty_table(){
    return new ElementTable{{ELEMENT_C, 0}, {ELEMENT_C13, 0}, {ELEMENT_H, 0}, {ELEMENT_H2, 0}, {ELEMENT_N, 0}, {ELEMENT_N15, 0}, {ELEMENT_O, 0}, {ELEMENT_O17, 0}, {ELEMENT_O18, 0}, {ELEMENT_P, 0}, {ELEMENT_P32, 0}, {ELEMENT_S, 0}, {ELEMENT_S34, 0}, {ELEMENT_S33, 0}, {ELEMENT_F, 0}, {ELEMENT_Cl, 0}, {ELEMENT_Br, 0}, {ELEMENT_I, 0}, {ELEMENT_As, 0}};
}



vector<string>* goslin::split_string(string text, char separator, char _quote, bool with_empty){
    bool in_quote = false;
    vector<string> *tokens = new vector<string>();
    stringstream sb;
    char last_char = '\0';
    bool last_escaped_backslash = false;
    
    for (uint32_t i = 0; i < text.length(); ++i){
        char c = text[i];
        bool escaped_backslash = false;
        if (!in_quote){
            if (c == separator){
                string sb_string;
                sb_string = sb.str();
                if (sb_string.length() > 0 || with_empty) tokens->push_back(sb_string);
                sb.str("");
            }
            else{
                if (c == _quote) in_quote = !in_quote;
                sb << c;
            }
        }
        else {
            if (c == '\\' && last_char == '\\' && !last_escaped_backslash){
                escaped_backslash = true;
            }
            else if (c == _quote && !(last_char == '\\' && !last_escaped_backslash)){
                in_quote = !in_quote;
            }
            sb << c;
        }
        
        last_escaped_backslash = escaped_backslash;
        last_char = c;
    }
    
    string sb_string;
    sb_string = sb.str();
    
    if (sb_string.length() > 0 || (last_char == ',' && with_empty)){
        tokens->push_back(sb_string);
    }
    if (in_quote){
        delete tokens;
        throw RuntimeException("Error: corrupted token in grammar");
    }
    return tokens;
}


string goslin::replace_all(std::string str, const std::string& from, const std::string& to) {
    int start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != (int)std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}


