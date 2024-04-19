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



template <class T> 
BaseParserEventHandler<T>::BaseParserEventHandler(){
    registered_events = new map<string, function<void(TreeNode *)>>();
    debug = "";
    error_message = "";
    word_in_grammar = false;
}

template <class T> 
BaseParserEventHandler<T>::~BaseParserEventHandler(){
    delete registered_events;
}

// checking if all registered events are reasonable and orrur as rules in the grammar
template <class T> 
void BaseParserEventHandler<T>::sanity_check(){
    
    for (auto event_name : *registered_events){
        
        
        if (!endswith(event_name.first, "_pre_event") && !endswith(event_name.first, "_post_event")){
            throw RuntimeException("Parser event handler error: event '" + event_name.first + "' does not contain the suffix '_pre_event' or '_post_event'");
        }
        
        string rule_name = event_name.first;
        rule_name = replace_all(rule_name, "_pre_event", "");
        rule_name = replace_all(rule_name, "_post_event", "");
        if (rule_names.find(rule_name) == rule_names.end()){
            throw RuntimeException("Parser event handler error: rule '" + rule_name + "' in event '" + event_name.first + "' is not present in the grammar");
        }
    }
}


template <class T> 
void BaseParserEventHandler<T>::handle_event(string event_name, TreeNode *node){
    if (debug == "full"){
        string reg_event = contains_val_p(registered_events, event_name) ? "*" : "";
        cout << event_name << reg_event << ": \"" << node->get_text() << "\"" << endl;
    }
        
    if (contains_val_p(registered_events, event_name)){
        if (debug != "" && debug != "full"){
            cout << event_name << ": \"" << node->get_text() << "\"" << endl;
        }
        registered_events->at(event_name)(node);
    }
}
