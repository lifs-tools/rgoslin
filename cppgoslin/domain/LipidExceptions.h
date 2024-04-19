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


#ifndef LIPID_EXCEPTIONS_H
#define LIPID_EXCEPTIONS_H

#include <string>

using namespace std;

class LipidException : public std::exception {
public:
    string message;
    LipidException(string _message){
        message = _message;
    }
    
    const char * what() const throw(){
        return message.c_str();
    }
};


class IllegalArgumentException : public LipidException {
public:
    IllegalArgumentException(string message) : LipidException("IllegalArgumentException: " + message){
        
    }
};


class ConstraintViolationException : public LipidException {
public:
    ConstraintViolationException(string message) : LipidException("ConstraintViolationException: " + message){
        
    }
};


class RuntimeException : public LipidException {
public:
    RuntimeException(string message) : LipidException("RuntimeException: " + message){
        
    }
};


class UnsupportedLipidException : public LipidException {
public:
    UnsupportedLipidException(string message) : LipidException("UnsupportedLipidException: " + message){
        
    }
};


class LipidParsingException : public LipidException {
public:
    LipidParsingException(string message) : LipidException("LipidParsingException: " + message){
        
    }
};


#endif /* LIPID_EXCEPTIONS_H */
