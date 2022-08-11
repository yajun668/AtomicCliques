#include <stdio.h>
#include "common.h"
#include <algorithm>
//This file is for common functions
//function: convert long to string
string itos (long a)
{
    ostringstream str;
    str<<a;
    return str.str();
}

//function: split Str
void SplitStr(vector<string> & splits, const string& str, const string & deli)
{	// Note: splits are never cleared inside!!!
    size_t start, end;

    start = str.find_first_not_of(deli);
    while (start != string::npos)
    {
        end = str.find_first_of(deli, start + 1);

        if (end != string::npos)
        {
            splits.push_back(str.substr(start, end - start));
            start = str.find_first_not_of(deli, end + 1);
        }
        else
        {
            splits.push_back(str.substr(start));
            break;
        }
    }
}


// find common vertices between two sorted subsets
vector<long> findCommonV(vector<long> &S1, vector<long> &S2)
{
    /* Must sort 2 subsets before call this function */
    sort(S1.begin(),S1.end());
    sort(S2.begin(),S2.end());
    vector<long> commonV;
    if (S1.empty() || S2.empty()) {
        return commonV;
    }

    long pos1=0;
    long pos2 =0;
    long p1,p2;
    long end1 = S1.size(), end2 = S2.size();
    while(pos1<end1 && pos2<end2)
    {
        p1 = S1[pos1];
        p2 = S2[pos2];
        if (p1 == p2) {
            commonV.push_back(p1);
            pos1++;
            pos2++;
        }
        else if (p1>p2)
            pos2++;
        else
            pos1++;
    }
    return commonV;
}
