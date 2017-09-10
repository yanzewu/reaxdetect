#pragma once

#include <vector>
#include <map>
#include <string>
#include "stringconvert.h"

using namespace std;

//unified interface of configuration
class ConfigReader {

public:
	ConfigReader() : split_char(' '){}
	explicit ConfigReader(const char split_char_) : split_char(split_char_){}
	string& operator[](const string& itemName) {
		return _items[itemName];
	}

	void read_data(const string& path) {
		fstream file(path, ios_base::in);
		if (!file.is_open()) {
			file.open(path, ios_base::out);
			for (const auto& item : _items) {
				file << item.first << split_char << item.second << endl;
			}
		}
		else {
			char buffer[256];
			while (file.getline(buffer, 256)) {
				auto split_result = split(buffer, -1, split_char);
				if (split_result.size() <= 1) {
					printf("Error: Invalid Argument for config input!\n");
					exit(1);
				}
				_items[split_result[0]] = split_result[1];
			}
		}
		file.close();
	}

	void fresh_data(const string& path)const {
		ofstream file(path, ios_base::out);
		for (const auto& item : _items) {
			file << item.first << split_char << item.second << endl;
		}
		file.close();
	}

	const string& get(const string& name, const string& d) {
		auto result = _items.find(name);
		if (result == _items.end())return d;
		else
			return result->second;
	}
	map<string, string> _items;
	const char split_char;
};