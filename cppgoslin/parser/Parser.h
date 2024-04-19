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


#ifndef PARSER_H
#define PARSER_H


#include "cppgoslin/parser/BaseParserEventHandler.h"
#include "cppgoslin/domain/StringFunctions.h"
#include "cppgoslin/parser/SumFormulaParserEventHandler.h"
#include "cppgoslin/parser/ParserClasses.h"
#include "cppgoslin/domain/Element.h"
#include "cppgoslin/parser/KnownGrammars.h"
#include <cstdint>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <math.h>

enum Content {NoContext, InLineComment, InLongComment, InQuote};
enum MatchWords {NoMatch, LineCommentStart, LineCommentEnd, LongCommentStart, LongCommentEnd, Quote};

using namespace std;

        

template <class T>
class Parser {
public:
    static const uint32_t SHIFT;
    static const uint64_t MASK;
    static const char RULE_ASSIGNMENT;
    static const char RULE_SEPARATOR;
    static const char RULE_TERMINAL;
    static const char EOF_SIGN;
    static const uint64_t EOF_RULE;
    static const uint64_t START_RULE;
    static const string EOF_RULE_NAME;
    
    
    uint64_t next_free_rule_index;
    map<char, set<uint64_t>> TtoNT;
    map<char, uint64_t> originalTtoNT;
    map<uint64_t, set<uint64_t>> NTtoNT;
    map<uint64_t, string> NTtoRule;
    map<uint64_t, vector<uint64_t>*> substitution;
    vector<Bitfield*> right_pair;
    int avg_pair;
    char quote;
    BaseParserEventHandler<T> *parser_event_handler;
    string grammar_name;
    bool used_eof;
    
    
    Parser(BaseParserEventHandler<T> *_parserEventHandler, string grammar_filename, char _quote = DEFAULT_QUOTE);
    Parser(BaseParserEventHandler<T> *_parserEventHandler, GrammarString grammar_string, char _quote = DEFAULT_QUOTE);
    void read_grammar(string grammar);
    virtual ~Parser();
    uint64_t get_next_free_rule_index();
    vector<string>* extract_text_based_rules(string grammar_filename, char _quote = DEFAULT_QUOTE);
    vector<uint64_t>* top_nodes(uint64_t rule_index);
    static uint64_t compute_rule_key(uint64_t rule_index_1, uint64_t rule_index_2);
    bool is_terminal(string product_token, char _quote);
    static string de_escape(string text, char _quote);
    uint64_t add_terminal(string text);
    vector<uint64_t>* collect_one_backwards(uint64_t rule_index);
    vector< vector<uint64_t>* >* collect_backwards(uint64_t child_rule_index, unsigned parent_rule_index);
    vector< vector<uint64_t>* >* collect_backwards(uint64_t child_rule_index, unsigned parent_rule_index, set<uint64_t>* visited, vector<uint64_t>* path, vector< vector<uint64_t>* >* collection);
    void fill_tree(TreeNode *node, DPNode *dp_node);
    string get_error_message();
    
    void raise_events(TreeNode *node);
    void raise_events_parallel(TreeNode *node, BaseParserEventHandler<T>*);
    
    virtual T parse(string text_to_parse, bool throw_error = true);
    virtual T parse_parallel(string text_to_parse, bool throw_error = true, BaseParserEventHandler<T>* bpeh = 0);
    
    TreeNode* parse_regular(string text_to_parse, BaseParserEventHandler<T>* bpeh = 0);
};



#include "cppgoslin/parser/Parser_impl.h"
#endif /* PARSER_H */


