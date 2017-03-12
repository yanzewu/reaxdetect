#include "reaxreader.h"

int ReaxReader::totmpscore(int weight, int nhcon, int terminalh) {
	return (weight * maxnhcon + nhcon) * maxhnum + terminalh;
}
int ReaxReader::tmpscore2score(int tmpscore) {
	return tmpscore % maxhnum + ((tmpscore / maxhnum) % maxnhcon * maxatmwht + (tmpscore / maxhnum) / maxnhcon) * maxhnum;
}
bool ReaxReader::isterminalhydrogen(int tmpscore) {
	return tmpscore == (maxnhcon + 1) * maxhnum;
}

void ReaxReader::MatMolecule::push_back(int index) {
	atoms.push_back(index);
}

bool ReaxReader::MatMolecule::operator==(const MatMolecule& othermlc)const {
	auto atom1 = atoms.begin();
	
	auto atom2 = othermlc.atoms.begin();
	while (atom1 != atoms.end() && atom2 != othermlc.atoms.end()) {
		if ((*pscores)[*atom1] != (*othermlc.pscores)[*atom2])return false;
		if ((pbondmatrix->begin() + *atom1)->size() != (othermlc.pbondmatrix->begin() + *atom2)->size())return false;
		atom1++; atom2++;
	}
	return true;

}
smiles ReaxReader::MatMolecule::to_smiles()const {
	//O(N)
	smiles sms;
	tsize_t i = 0;
	Array atomrpos(pbondmatrix->size(), -2);
	for (const auto& atom : atoms) {
		if (isterminalhydrogen((*pscores)[atom])) {
			atomrpos[atom] = -1;
			continue;//no terminal hydrogen
		}
		else if ((*pscores)[atom] == maxnhcon * maxhnum + 1) {	// (X-)H-H
			if (atomrpos[(*pbondmatrix)[atom].front()] == -2) {
				atomrpos[atom] = -1;
				continue;
			}
		}
		sms._push_back(tmpscore2score((*pscores)[atom]));
		atomrpos[atom] = i;
		i++;
	}
	//O(N*k)
	i = 0;
	for (const auto& atom : atoms) {
		if (atomrpos[atom] < 0)continue;
		for (const auto& child : (*pbondmatrix)[atom]) {
			if (atomrpos[child] >= 0)sms._push_child(i, atomrpos[child]);
		}
		i++;
	}
	sms.canonicalize();
	return sms;
}
