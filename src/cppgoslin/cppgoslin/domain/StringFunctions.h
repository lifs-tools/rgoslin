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


#ifndef STRINGFUNCTIONS_H
#define STRINGFUNCTIONS_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/Element.h"
#include <algorithm>
#include <cctype>
#include <cstdint>


#define contains_val(X, Y) ((X).find(Y) != (X).end())
#define contains_val_p(X, Y) ((X)->find(Y) != (X)->end())
#define uncontains_val(X, Y) ((X).find(Y) == (X).end())
#define uncontains_val_p(X, Y) ((X)->find(Y) == (X)->end())

static const char DEFAULT_QUOTE = '\'';

using namespace std;

namespace goslin {
    string replace_all(std::string str, const std::string& from, const std::string& to);
    string strip(string s, char c);
    vector<string>* split_string(string text, char separator, char _quote = DEFAULT_QUOTE, bool with_empty = false);
    string to_lower(string s);
    string to_upper(string s);
    bool endswith(const string &main_str, const string &to_match);
    double get_mass(ElementTable *elements);
}

#endif /* STRINGFUNCTIONS_H */
