#include "reaxreader.h"
#include "listhandler.h"

void find_neighbour(int index, int root, ReaxReader::Matrix* bondmat, ReaxReader::Mark skip_atoms, ReaxReader::Mark mark) {
	mark[index] = true;
	for (auto child = bondmat->at(index).begin(); child != bondmat->at(index).end(); child++) {
		if (mark[*child])continue;
		if (!skip_atoms[*child]) {
			mark[*child] = true;
			bondmat->at(root).push_back(*child);
		}
		else {
			find_neighbour(*child, root, bondmat, skip_atoms, mark);
		}
	}
}

void ReaxReader::_skip_bond(ReaxReader::Matrix* bondmat) {

	Mark mark = new bool[bondmat->size()]{ false };

	for (int i = 1; i < bondmat->size(); i++) {
		if (skip_atoms[i]) continue;
		memset(mark, 0, sizeof(bool)*bondmat->size());
		mark[i] = true;
		for (auto child = bondmat->at(i).begin(); child != bondmat->at(i).end(); child++) {
			if (!skip_atoms[*child])mark[*child] = true;
		}

		for (auto child = bondmat->at(i).begin(); child != bondmat->at(i).end();) {
			if (skip_atoms[*child]) {
				find_neighbour(*child, i, bondmat, skip_atoms, mark);
				child = bondmat->at(i).erase(child);
			}
			else {
				child++;
			}
		}
	}
	for (int i = 1; i < bondmat->size(); i++) {
		if (skip_atoms[i]) {
			bondmat->at(i).clear();
		}
	}
	delete[] mark;
}
