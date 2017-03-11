#pragma once

#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

string str(double num);

int stob(const string&);

int char2int(const char c);
double string2double(const string& str);

vector<string> split(const string& str, int num = -1, char symbol = ' ');

string join(const vector<string>&, const string& symbol = ",");
string join(const vector<int>&, const string& symbol = ",");
string join(const vector<double>&, const string& symbol = ",");
