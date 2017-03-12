#include <string>
#include "reaxreader.h"

bool ReaxReader::Reaction::check_valid()
{
	sort(reagants.begin(), reagants.end(), less<int>());
	sort(products.begin(), products.end(), less<int>());

	return (reagants != products);
}
ReaxReader::Reaction ReaxReader::Reaction::operator-()const {
	return ReaxReader::Reaction(products, reagants);
}

bool ReaxReader::Reaction::operator==(const Reaction & reac) const
{
	return (reagants == reac.reagants) & (products == reac.products);
}
string ReaxReader::Reaction::to_string(const vector<string>& species)const {
	string out;
	vector<int>::const_iterator iter;
	for (iter = reagants.begin(); iter < reagants.end() - 1; iter++) {
		out.append(species[*iter] + "+");
	}
	out.append(species[*iter]+ "=");
	for (iter = products.begin(); iter < products.end() - 1; iter++) {
		out.append(species[*iter] + "+");
	}
	out.append(species[*iter]);
	return out;
}
