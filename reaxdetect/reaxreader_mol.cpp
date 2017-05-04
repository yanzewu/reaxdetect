#include "reaxreader.h"
#include "listhandler.h"

int ReaxReader::totmpscore(int weight, int nhcon, int bond, int terminalh) {
	return ((weight * maxnhcon + nhcon) * maxbond + bond) * maxhnum + terminalh;
}
int ReaxReader::tmpscore2score(int tmpscore) {
	return tmpscore % (maxhnum * maxbond)+ 
		((tmpscore / (maxhnum * maxbond)) % maxnhcon * maxatmwht + (tmpscore / (maxhnum * maxbond)) / maxnhcon) * maxbond * maxhnum;
}
bool ReaxReader::isterminalhydrogen(int tmpscore) {
	return tmpscore == ((maxnhcon + 1) * maxbond + 1) * maxhnum;
}

void ReaxReader::MatMolecule::push_back(int index) {
	atoms.push_back(index);
}

bool ReaxReader::MatMolecule::operator==(const MatMolecule& othermlc)const {
	if (atoms.size() != othermlc.atoms.size())return false;
	
	auto atom1 = atoms.begin();
	
	auto atom2 = othermlc.atoms.begin();
	while (atom1 != atoms.end() && atom2 != othermlc.atoms.end()) {
		if ((*pscores)[*atom1] != (*othermlc.pscores)[*atom2])return false;
		if ((pbondmatrix->begin() + *atom1)->size() != (othermlc.pbondmatrix->begin() + *atom2)->size())return false;
		atom1++; atom2++;
	}
	return true;

}
bool ReaxReader::MatMolecule::equals_to(MatMolecule & othermlc)
{
	return contain_equal(atoms, othermlc.atoms, less<int>(), equal_to<int>());
}
smiles ReaxReader::MatMolecule::to_smiles()const {
	//O(N)
	smiles sms;
	tsize_t i = 0;
	Array atomrpos(pbondmatrix->size(), -2);
	for (const auto& atom : atoms) {
		if (isterminalhydrogen((*pscores)[atom]) && (*pbondmatrix)[atom].front().second != 0) {
			atomrpos[atom] = -1;
			continue;//no terminal hydrogen
		}
		else if ((*pscores)[atom] == (maxnhcon * maxbond + 1) * maxhnum + 1) {	// (X-)H-H
			if (atomrpos[(*pbondmatrix)[atom].front().first] == -2) {
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
			if (atomrpos[child.first] >= 0)sms._push_child(i, atomrpos[child.first], child.second);
		}
		i++;
	}
	sms.canonicalize();
	return sms;
}
