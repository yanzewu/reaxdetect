#pragma once

#include <vector>
#include <map>
#include <string>
#include "stringconvert.h"

using namespace std;

//unified interface of configuration
class ConfigReader {

public:
	string& operator[](string itemName) {
		return _items[itemName];
	}

	void read_data(const string& path) {
		fstream file(path, ios_base::in);
		if (!file.is_open()) {
			file.open(path, ios_base::out);
			for (const auto& item : _items) {
				file << item.first << " " << item.second << endl;
			}
		}
		else {
			char buffer[256];
			while (file.getline(buffer, 256)) {
				auto split_result = split(buffer);
				if (split_result.size() <= 1)runtime_error("Config: Invalid Argument!");
				_items[split_result[0]] = split_result[1];
			}
		}
		file.close();
	}

	string get(const string& name, const string& d) {
		auto result = _items.find(name);
		if (result == _items.end())return d;
		else
			return result->second;
	}
private:
	map<string, string> _items;
};