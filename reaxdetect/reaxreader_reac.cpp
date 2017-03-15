#include <string>
#include "reaxreader.h"
#include "listhandler.h"

bool ReaxReader::Reaction::check_valid()
{
	return !contain_equal(reagants, products, less<int>(), equal_to<int>());
}
ReaxReader::Reaction ReaxReader::Reaction::operator-()const {
	return ReaxReader::Reaction(products, reagants);
}

bool ReaxReader::Reaction::operator==(const Reaction & reaction) const
{
	return (reagants == reaction.reagants) & (products == reaction.products);
}
bool ReaxReader::Reaction::operator!=(const Reaction & reaction) const
{
	return (reagants != reaction.reagants) || (products != reaction.products);
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
