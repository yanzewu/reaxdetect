#pragma once

#ifndef PATH_H
#define PATH_H

#include <string>
#include <vector>

using namespace std;

namespace path {
	bool exists(const string&);
	string get_path(const string&);
	vector<string> split(const string&);
}

#endif