#pragma once
#include "includes.h"
#include "typedef.h"

/* Contains templates for some quick operations */

//find if a data is in a vector, if not then push_back. also adds frequency
template<typename type, typename container>
size_t pushnew(vector<type>& datas, const type& data, container& freq) {
	size_t pos = find(datas.begin(), datas.end(), data) - datas.begin();
	if (pos == datas.size()) {
		datas.push_back(data);
		freq.push_back(1);
		return pos;
	}
	else {
		freq[pos]++;
		return pos;
	}
}

//find if a data is in a vector, if not then push_back. also adds frequency and records the index.
template<typename type>
size_t pushnew2index(unsigned int index, const vector<type>& datas, vector<unsigned int>& indexes, vector<int>& freq) {

	for (unsigned int i = 0; i < indexes.size(); i++) {
		if (datas[indexes[i]] == datas[index]) {
			freq[i]++; 
			return i;
		}
	}
	indexes.push_back(index);
	freq.push_back(1);
	return indexes.size() - 1;
}

// will change content by sorting.
template<typename container, typename order_sort, typename order_equal>
bool contain_equal(container& c1, container& c2, order_sort pred_sort, order_equal pred_equal) {
	sort(c1.begin(), c1.end(), pred_sort);
	sort(c2.begin(), c2.end(), pred_sort);
	return (c1.size() == c2.size()) && (equal(c1.begin(), c1.end(), c2.begin(), pred_equal));
}

// will change content by sorting.
template<typename type, typename order_sort, typename order_equal>
bool contain_equal(list<type>& c1, list<type>& c2, order_sort pred_sort, order_equal pred_equal) {
	c1.sort(pred_sort);
	c2.sort(pred_sort);
	return (c1.size() == c2.size()) && (equal(c1.begin(), c1.end(), c2.begin(), pred_equal));
}

//sort a index list by score. Pred means the order of head to tail.
template<typename type, typename order>
void sortby(list<type>& source, const Array& score, order pred) {
	if (source.empty())return;

	auto _pred = [score, pred](auto x, auto y) {
		return pred(score[x], score[y]);
	};
	source.sort(_pred);
}

//sort a index list by score. Pred means the order of head to tail.
template<typename container, typename order>
void sortby(container& source, const Array& score, order pred) {
	if (source.empty())return;

	auto _pred = [score, pred](auto x, auto y) {
		return pred(score[x], score[y]);
	};
	sort(source.begin(), source.end(), _pred);
}

//rank by a exist rank
template<typename iter>
void sortwithrank(iter begin, iter end, Array& rank) {
	if (end == begin || end == begin + 1)return;
	unsigned int i = 0;
	for (iter s = begin; s != end; s++) {
		if (rank[i] != i) {
			swap(*s, *(begin + rank[i]));
			int t = rank[rank[i]];
			rank[rank[i]] = rank[i];
			rank[i] = t;
		}
		i++;
	}
}

//sort a list with another list accompanied. Pred means the order of head to tail.
template<typename type1, typename type2, typename order>
void sortwith(vector<type1>& score, vector<type2>& assist, const order pred) {
	for (size_t i = 0; i < score.size() - 1; i++) {
		for (size_t j = 0; j < score.size() - 1 - i; j++) {
			if (pred(score[j + 1], score[j])) {
				swap(score[j], score[j + 1]);
				swap(assist[j], assist[j + 1]);
			}
		}
	}

}
