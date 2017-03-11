#pragma once

#include <utility>
#include <vector>
#include "typedef.h"

using namespace std;

//sparce vector
template<typename type>
class svector {
public:

	svector() : svector(0) {
	}
	explicit svector(type empty) : _size(0), empty(empty) {
	}
	size_t size()const {
		return _size;
	}
	void resize(size_t newsize) {
		_size = newsize;
	}

	type& operator[](const tsize_t pos) {
		auto findpos = find(poses.begin(), poses.end(), pos);
		if (findpos == poses.end()) {
			push_back(empty, pos);
			return datas.back();
		}
		else {
			return datas[findpos - poses.begin()];
		}
	}
	const type& operator[](tsize_t pos)const {
		auto findpos = find(poses.begin(), poses.end(), pos);
		if (findpos == poses.end())return empty;
		else {
			return datas[findpos - poses.begin()];
		}
	}
	void push_back(const type& data) {
		if (data == empty) {
			_size++;
		}
		else {
			poses.push_back(++_size);
			datas.push_back(data);
		}
	}
	void push_back(const type& data, const tsize_t pos) {
		poses.push_back(pos);
		datas.push_back(data);
		_size = max(_size, pos + 1);
	}
	vector<type> to_array()const {
		vector<type> out_array(_size, empty);
		for (size_t i = 0; i < poses.size(); i++) {
			out_array[poses[i]] = datas[i];
		}
		return out_array;
	}

	vector<tsize_t> poses;
	vector<type> datas;
protected:
	tsize_t _size;
	type empty;
};


template<typename type>
size_t pushnewpm(vector<type>& datas, const type& data, svector<diint>& freq) {

	int cnt = 0;
	for (auto datas_ = datas.begin(); datas_ != datas.end(); datas_++) {
		if (*datas_ == data) {
			freq[datas_ - datas.begin()].first++;
			return datas_ - datas.begin();
		}
		else if (*datas_ == -data) {
			freq[datas_ - datas.begin()].second++;
			return datas_ - datas.begin();
		}
	}
	datas.push_back(data);
	freq.push_back(diint(1, 0));
	return datas.size();
}

