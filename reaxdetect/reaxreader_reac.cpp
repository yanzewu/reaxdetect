#include <string>
#include <algorithm>

#include "reaxreader.h"
#include "util/algorithmutil.h"
#include "util/strutil.h"

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

	vector<string> reactants_str, products_str;
	for (const auto& reactant : reagants){
		reactants_str.push_back(species[reactant]);
	}
	sort(reactants_str.begin(), reactants_str.end());
	for (const auto& product : products) {
		products_str.push_back(species[product]);
	}
	sort(products_str.begin(), products_str.end());
	return join(reactants_str, "+") + "=" + join(products_str, "+");
}
