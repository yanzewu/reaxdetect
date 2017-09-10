#pragma once

#ifndef VECUTIL_H
#define VECUTIL_H

#include <vector>

using namespace std;

namespace veccal {
	template<typename _Type>
	vector<_Type> operator+(const vector<_Type>& a, const vector<_Type>& b) {
		vector<_Type> c(a);
		auto c_ = c.begin();
		auto b_ = b.begin();
		for (; b_ < b.end(); b_++, c_++) {
			*c_ += *b_;
		}
		return c;
	}

	template<typename _Type>
	void operator+=(vector<_Type>& a, const vector<_Type>& b) {
		auto a_ = a.begin(); auto b_ = b.begin();
		for (; b_ < b.end(); b_++, a_++) {
			*a_ += *b_;
		}
	}

	template<typename _Type>
	void operator*=(vector<_Type>& a, double b) {
		for (auto& a_ : a)a_ *= b;
	}

	template<typename _Type>
	vector<_Type> operator*(const vector<_Type>& a, double b) {
		vector<_Type> dst(a);
		dst *= b;
		return dst;
	}

	template<typename _Type_src, typename _Type_dst>
	void mul_into(const vector<_Type_src>& src, double d, vector<_Type_dst>& dst) {
		dst.resize(src.size());
		
		auto src_ = src.begin();
		auto dst_ = dst.begin();
		for (; src_ < src.end(); src_++, dst_++) {
			*dst_ = (*src_) * d;
		}
	}

	template<typename _Type_src, typename _Type_dst, typename _Type_d>
	void multiply(vector<_Type_src>& src, const vector<_Type_d>& d, const vector<_Type_dst>& dst) {
		dst.resize(src.size());

		auto src_ = src.begin();
		auto d_ = d.begin();
		auto dst_ = dst.begin();

		for (; src_ < src.end(); src_++, d_++, dst_++) {
			*dst_ = (*src_) * (*d_);
		}
	}

	template<typename _Type_src, typename _Type_dst, typename _Type_d>
	void divide(vector<_Type_src>& src, const vector<_Type_d>& d, vector<_Type_dst>& dst) {
		dst.resize(src.size());

		auto src_ = src.begin();
		auto d_ = d.begin();
		auto dst_ = dst.begin();

		for (; src_ < src.end(); src_++, d_++, dst_++) {
			*dst_ = (*src_) / *d_;
		}
	}

	//find the minimum one in an array.
	template<typename iter>
	iter mininum(const iter first, const iter last) {
		iter m = first;
		for (iter p = first; p != last; p++) {
			if (*p < *m)m = p;
		}
		return m;
	}

	//find the maximum one in an array
	template<typename iter>
	iter findmax(const iter first, const iter last) {
		iter m = first;
		for (iter p = first; p != last; p++) {
			if (*p > *m)m = p;
		}
		return m;
	}

}

template<typename _Type>
vector<_Type> net(const vector<pair<_Type, _Type> >& src) {
	vector<_Type> dst(src.size(), 0);
	auto src_ = src.begin();
	auto dst_ = dst.begin();
	for (; src_ < src.end(); src_++, dst_++) {
		*dst_ = src_->first - src_->second;
	}
	return dst;
}

template<typename _Type_src, typename _Type_dst>
void veccpy(const vector<_Type_src>& src, vector<_Type_dst>& dst) {
	dst.resize(src.size());
	auto src_ = src.begin();
	auto dst_ = dst.begin();
	for (; src_ < src.end(); src_++, dst_++) {
		*dst_ = (_Type_dst)*src_;
	}
}

template<typename _Type>
void vwrite(ofstream& outfile, const vector<_Type>& v) {
	unsigned int size = v.size();
	outfile.write((char*)&size, sizeof(size));
	if (size == 0)return;
	outfile.write((char*)&v[0], size * sizeof(_Type));
}

template<typename _Type>
void vread(ifstream& infile, vector<_Type>& v) {
	unsigned int size = 0;
	infile.read((char*)&size, sizeof(size));
	if (size == 0)return;
	v.resize(size);
	infile.read((char*)&v[0], size * sizeof(_Type));
}

#endif