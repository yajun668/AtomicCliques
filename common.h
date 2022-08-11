#ifndef common_h
#define common_h
#endif /* common_h */
#include <vector>
#include <string>
#include <iostream>
#include<fstream>
#include <sstream>
using namespace std;
//This file is for common functions
string itos (long a);
void SplitStr(vector<string> & splits, const string & str, const string& deli = " "); //function: split Str
vector<long> findCommonV(vector<long> &S1, vector<long> &S2);





