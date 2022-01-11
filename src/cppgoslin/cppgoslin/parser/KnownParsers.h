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


#ifndef KNOWN_PARSERS_H
#define KNOWN_PARSERS_H


#include "cppgoslin/parser/FattyAcidParserEventHandler.h"
#include "cppgoslin/parser/ShorthandParserEventHandler.h"
#include "cppgoslin/parser/GoslinParserEventHandler.h"
#include "cppgoslin/parser/LipidMapsParserEventHandler.h"
#include "cppgoslin/parser/SwissLipidsParserEventHandler.h"
#include "cppgoslin/parser/HmdbParserEventHandler.h"
#include "cppgoslin/parser/KnownGrammars.h"

class FattyAcidParser : public Parser<LipidAdduct*> {
public:
    FattyAcidParser();
    ~FattyAcidParser();
    LipidAdduct* parse(string lipid_name, bool throw_error = true);
    LipidAdduct* parse_parallel(string, bool throw_error = true, BaseParserEventHandler<LipidAdduct*>* bpeh = 0);
};


class ShorthandParser : public Parser<LipidAdduct*> {
public:
    ShorthandParser();
    ~ShorthandParser();
    LipidAdduct* parse_parallel(string, bool throw_error = true, BaseParserEventHandler<LipidAdduct*>* bpeh = 0);
};


class GoslinParser : public Parser<LipidAdduct*> {
public:
    GoslinParser();
    ~GoslinParser();
    LipidAdduct* parse_parallel(string, bool throw_error = true, BaseParserEventHandler<LipidAdduct*>* bpeh = 0);
};


class LipidMapsParser : public Parser<LipidAdduct*> {
public:
    LipidMapsParser();
    ~LipidMapsParser();
    LipidAdduct* parse_parallel(string, bool throw_error = true, BaseParserEventHandler<LipidAdduct*>* bpeh = 0);
};


class SwissLipidsParser : public Parser<LipidAdduct*> {
public:
    SwissLipidsParser();
    ~SwissLipidsParser();
    LipidAdduct* parse_parallel(string, bool throw_error = true, BaseParserEventHandler<LipidAdduct*>* bpeh = 0);
};


class HmdbParser : public Parser<LipidAdduct*> {
public:
    HmdbParser();
    ~HmdbParser();
    LipidAdduct* parse_parallel(string, bool throw_error = true, BaseParserEventHandler<LipidAdduct*>* bpeh = 0);
};


class LipidParser {
public:
    vector<Parser<LipidAdduct*>*> parser_list;
    Parser<LipidAdduct*>* lastSuccessfulParser;
    
    LipidParser();
    ~LipidParser();
    LipidAdduct* parse(string lipid_name);
    LipidAdduct* parse_parallel(string lipid_name);
};      

#endif /* KNOWN_PARSERS_H */
