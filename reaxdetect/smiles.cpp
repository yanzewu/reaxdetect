//a newer version of molecule
#include"smiles.h"
#include<stack>
#include"stringconvert.h"
#include"listhandler.h"
#include"elementname.h"

#define lettermode	1
#define hmode		2

const int primes[35] = {
	2,3,5,7,11, 13,17,19,23,29, 31,33,37,43,47, 53,59,61,67,71,
	73,79,83,89,97, 101,103,107,109,113, 127,131,137,139,149
};

bool smiles::operator==(const smiles& sms)const {
	if (atoms.size() != sms.atoms.size())return false;
	bool* mark = new bool[atoms.size()]{ false };
	return compare(sms, 0, 0, mark);
}
bool smiles::compare(const smiles& sms, int myindex, int hisindex, bool* mark)const {
	mark[myindex] = true;
	if (atoms[myindex].score != sms.atoms[hisindex].score)return false;
	if (atoms[myindex].children.size() != sms.atoms[hisindex].children.size())return false;
	for (auto mychild = atoms[myindex].children.begin(), hischild = sms.atoms[hisindex].children.begin();
		mychild != atoms[myindex].children.end(); mychild++, hischild++) {
		if (mark[*mychild])continue;
		if (!compare(sms, *mychild, *hischild, mark))return false;
	}
	return true;
}
void smiles::localconvert()
{
	//1.copy data (score and childtmp is modifyable)
	vector<list<int> > _children_tmp(atoms.size());//neighbour
	vector<int> _score(atoms.size()), rank(atoms.size());//score is the unstable score; no matter what happens, score is host. rank is client.
	for (unsigned int i = 0; i < atoms.size(); i++) {
		_children_tmp[i] = atoms[i].children;
		_score[i] = atoms[i].score;
	}

	//2.get symmmetry classes (O(N^3))
	while (rankby(_score, rank)) {//if rank is changed
		//O(N*k)
		for (unsigned int i = 0; i < _score.size(); i++) {
			renew_score(i, rank, _children_tmp[i], _score);
		}
	}
	
	//take atom weight into account...
	for (unsigned int i = 0; i < _score.size(); i++) {
		_score[i] = rank[i] * maxatmwht * maxhnum * maxnhcon + atoms[i].score;
	}
	rankby(_score, rank);

	//3.canonical labeling
	//find the start
	_score = rank;//here score is served as the symmetry class. (host) rank is served as label. (client)

	int lowest = 0, crtlabel = 1;
	for (const auto& s : _score) {
		if (s > 0)lowest++;
		else break;
	}
	for (auto& r : rank)
		r = -1;
	rank[lowest] = 0;//clear the label vector.

	//tie break algorithm from open babel
	label(lowest, crtlabel, _children_tmp, _score, rank);

	//4.O( Nlog(N) ) rank atoms
	auto rank_cpy = rank;
	sortwith(rank_cpy, atoms, less<int>());//put the order of atoms

	//O(N*m) (replace children and rank children, as at this time children is the old index)
	for (auto& atom : atoms) {
		for (auto& child : atom.children)
			child = rank[child];
		atom.children.sort(greater<int>());
	}
}
bool smiles::rankby(vector<int>& score, vector<int>& rank) {
	//a more complex algorithm...
	//first we rank the array, then...

	vector<int> sort_result(score.size());
	for (unsigned int i = 0; i < score.size(); i++)
		sort_result[i] = i;
	sortby(sort_result, score, less<int>());

	//compare and delete equal rank
	//backup old rank
	vector<int> _newrank(rank.size());

	int highest_in_old = 0;
	for (const auto& r : rank)
		if (r > highest_in_old)highest_in_old = r;

	_newrank[sort_result[0]] = 0;
	int rank_cnt = 0;
	for (unsigned int i = 1; i < rank.size(); i++) {
		if (score[sort_result[i]] > score[sort_result[i - 1]])rank_cnt++;
		_newrank[sort_result[i]] = rank_cnt;
	}

	if (rank_cnt > highest_in_old) {
		//copy rank
		rank = _newrank; return true;
	}
	else return false;
}
bool smiles::renew_score(const int index, const vector<int>& rank, const list<int>& _children_tmp, vector<int>& score) {
	int oldscore = score[index];
	score[index] = 1;
	for (auto child : _children_tmp) {
		score[index] *= primes[rank[child]];
	}
	return score[index] >= oldscore;
}
void smiles::label(int index, int& crtlabel, vector<list<int> >& children, const vector<int>& symclass, vector<int>& labels) {
	//remove already labeled
	for (auto child = children[index].begin(); child != children[index].end();) {
		if (labels[*child] >= 0)child = children[index].erase(child);
		else child++;
	}
	sortby(children[index], symclass, greater<int>());

	//then, label children
	for (const auto& child : children[index]) {
		labels[child] = crtlabel;
		crtlabel++;
	}

	//finally, label others
	for (const auto& child : children[index])
		label(child, crtlabel, children, symclass, labels);
}

int nameL2weight(const string& name) {
	size_t result = find(WeightNameL, WeightNameL + 32, name) - WeightNameL;
	if (result < 32)return result;
	else return 32;
}

int smiles::score2weight(int score) {
	return (score / maxhnum) % maxatmwht;
}
int smiles::score2hydro(int score) {
	return score % maxhnum;
}


int smiles::has_element(int weight)const {
	int out = 0;
	for (const auto& a : atoms) {
		if (score2weight(a.score) == weight)out++;
	}
	return out;
}

void smiles::compile()
{
	//scan loop&delete parent
	bool* mark = new bool[atoms.size()]{ 0 };

	//	list<bond> loops;
	scanloop(0, -1, mark);
	vector<string> loopmark(atoms.size());

	int crtloopid = 1;
	for (unsigned int i = 0; i < looproots.size(); i++)
		loopmark[looproots[i]] += to_string(i + 1);

	//compile name
	symbol = "";
	compilepart(0, loopmark);
	if (symbol == "hh-hh")symbol = "h2";
	delete[] mark;
}
void smiles::scanloop(int index, int parent, bool* mark)
{
	mark[index] = true;
	for (auto child = atoms[index].children.begin(); child != atoms[index].children.end();) {
		if (*child == parent) {
			child = atoms[index].children.erase(child);
		}
		else {
			if (mark[*child]) {
				atoms[*child].children.erase(
					find(atoms[*child].children.begin(), atoms[*child].children.end(), index));
				looproots.push_back(*child);
				*child = -(int)looproots.size();
			}
			else scanloop(*child, index, mark);
			child++;
		}
	}
}
void smiles::compilepart(int index, const vector<string>& loopmark) {
	//roots
	if (index < 0) {
		symbol.append(to_string(-index));
		return;
	}

	symbol.append(WeightNameL[score2weight(atoms[index].score)]);
	symbol.append(loopmark[index]);
	if (score2hydro(atoms[index].score) > 1)symbol.append("h" + to_string(score2hydro(atoms[index].score)));
	else if (score2hydro(atoms[index].score) == 1)symbol.append("h");

	for (auto child : atoms[index].children) {
		if (child != atoms[index].children.back()) {
			symbol.append("(");
			compilepart(child, loopmark);
			symbol.append(")");
		}
		else {
			symbol.append("-");
			compilepart(child, loopmark);
		}
	}
}

unsigned int smiles::toBin(int* buffer)const {
	int* p = buffer;
	*p = atoms.size(); p++;
	for (const auto& atom : atoms) {
		*p = atom.score; p++;
		*p = atom.children.size(); p++;
		for (const auto& child : atom.children) {
			*p = child; p++;
		}
	}
	*p = looproots.size(); p++;
	for (const auto& looproot : looproots) {
		*p = looproot; p++;
	}
	return (p - buffer) * sizeof(int);
}
void smiles::readBin(int* buffer) {
	int* p = buffer;
	atoms.resize(*p); p++;
	for (auto& atom : atoms) {
		atom.score = *p; p++;
		atom.children.resize(*p); p++;
		for (auto& child : atom.children) {
			child = *p; p++;
		}
	}
	looproots.resize(*p); p++;
	for (auto& looproot : looproots) {
		looproot = *p; p++;
	}

	vector<string> loopmark(atoms.size());

	int crtloopid = 1;
	for (unsigned int i = 0; i < looproots.size(); i++)
		loopmark[looproots[i]] += to_string(i + 1);

	compilepart(0, loopmark);
}

void smiles::readSmiles(const string& symbol) {
	//first read the atom, then scan the connectivity
	auto s = symbol.begin();

	struct bond { int id_1, id_2; };
	vector<bond> loops;
	stack<size_t> atomstack;
	string name_buffer;
	int prev_cnt = -1;
	char mode = lettermode;
	while (s != symbol.end()) {
		if (*s >= 65 && *s != 'h' && mode == lettermode)name_buffer.push_back(*s);
		else if (*s < 65 && *s >= 48) {
			//end letter
			if (mode == lettermode) {
				mode = 0;
				atomstack.push(atoms.size());
				atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
				name_buffer.clear();
				if (prev_cnt >= 0) {
					atoms.back().children.push_back(prev_cnt);
					atoms[prev_cnt].children.push_back(atoms.size() - 1);
				}
				goto pushloop;
			}
			//loops
			else if (mode == 0) {
			pushloop:
				if (loops.size() < char2int(*s)) {
					loops.resize(char2int(*s)); loops[char2int(*s)].id_1 = atoms.size() - 1;
				}
				else loops[char2int(*s)].id_2 = atoms.size() - 1;
			}
			//h
			else if (mode == hmode) {
				atoms.back().score += char2int(*s);
				mode = 0;
			}
		}
		else if (*s == 'h') {
			if (mode == lettermode) {
				atomstack.push(atoms.size());
				atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
				name_buffer.clear();
				if (prev_cnt >= 0) {
					atoms.back().children.push_back(prev_cnt);
					atoms[prev_cnt].children.push_back(atoms.size() - 1);
				}
			}
			mode = hmode;
		}
		else if (*s == '-' || *s == '(' || *s == ')') {
			if (mode == lettermode) {
				atomstack.push(atoms.size());
				atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
				name_buffer.clear();
				if (prev_cnt >= 0) {
					atoms.back().children.push_back(prev_cnt);
					atoms[prev_cnt].children.push_back(atoms.size() - 1);
				}
			}
			else if (mode == hmode) {
				atoms.back().score += 1; mode = 0;
			}
			prev_cnt = atomstack.top(); mode = lettermode;
			if (*s != '(')atomstack.pop();
		}
		s++;
	}
	if (mode == lettermode) {
		atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
		name_buffer.clear();
		if (prev_cnt >= 0) {
			atoms.back().children.push_back(prev_cnt);
			atoms[prev_cnt].children.push_back(atoms.size() - 1);
		}
	}
	else if (mode == hmode) {
		atoms.back().score += 1; mode = 0;
	}

	for (auto loop : loops) {
		atoms[loop.id_1].children.push_back(loop.id_2);
		atoms[loop.id_2].children.push_back(loop.id_1);
	}
	for (auto atom = atoms.begin(); atom != atoms.end(); atom++) {
		atom->score += maxhnum*maxatmwht*atom->children.size();
	}
	localconvert();
}