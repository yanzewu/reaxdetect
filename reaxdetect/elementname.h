#pragma once
#include"includes.h"

//used as element name storage

#define MAX_ELEMENT_SUPPORT	16
#define MAX_WEIGHT_SUPPORT	64

#define UNKNOWN_WEIGHT	"X"

const string ElementName[MAX_ELEMENT_SUPPORT] = {
	"h","He","Li","Be","B", "c","n","o","F","Ne", "Na","Mg","Al","Si","p", "s"
};
const string ElementNameL[MAX_ELEMENT_SUPPORT] = {
	"h","[He]","[Li]","[Be]","[B]", "c","n","o","[F]","[Ne]", "[Na]","[Mg]","[Al]","[Si]","p", "s"
};

const string WeightName[MAX_WEIGHT_SUPPORT] = {
	"", "h","D","T","He","X", "Li","Li","Be","X","X", "B","c","X","n","X", "o","X","X","F","Ne",
	"X","X","Na","Mg","X", "X","Al","Si","X","X", "p","s"
};

const string WeightNameL[MAX_WEIGHT_SUPPORT] = {
	"", "h","[D]","[T]","[He]",UNKNOWN_WEIGHT, "[Li]","[Li]","[Be]","X","X", "[B]","c","X","n","X", "o","X","X","F","Ne",
	"X","X","[Na]","[Mg]","X", "X","[Al]","[Si]","X","X", "p","s"
};
