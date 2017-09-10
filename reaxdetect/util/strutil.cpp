
#include "strutil.h"

vector<string> split(const string& str, int num, char symbol) {
	int cnt = 0;
	vector<string> out;
	int left = 0, right = 0;
	for (right = 0; right < str.size(); right++) {
		if (str[right] == symbol) {
			out.push_back(str.substr(left, right - left));
			left = right + 1;
			cnt++; if (cnt == num)break;
		}
	}
	out.push_back(str.substr(left, str.size() - left));
	return out;
}

string join(const vector<string>& src, const string& symbol) {
	if (src.empty()) {
		return "";
	}
	string dst = src[0];
	for (auto src_ = src.begin() + 1; src_ < src.end(); src_++) {
		dst.append(symbol + *src_);
	}
	return dst;
}
string join(const vector<int>& src, const string& symbol) {
	if (src.empty()) {
		return "";
	}
	string dst = str(src[0]);
	for (auto src_ = src.begin() + 1; src_ < src.end(); src_++) {
		dst.append(symbol + str(*src_));
	}
	return dst;
}
string join(const vector<double>& src, const string& symbol) {
	if (src.empty()) {
		return "";
	}
	string dst = str(src[0]);
	for (auto src_ = src.begin() + 1; src_ < src.end(); src_++) {
		dst.append(symbol + str(*src_));
	}
	return dst;
}


int char2int(const char c) {
	return c - 48;
}
double string2double(const string& str) {
	return strtod(str.c_str(), NULL);
}

int stob(const string& s) {
	if (s == "true") {
		return 1;
	}
	else if (s == "false") {
		return 0;
	}
	else {
		return -1;
	}
}

string str(double num) {
	char buffer[6];
	snprintf(buffer, sizeof(buffer), "%f", num);
	return string(buffer);
}
