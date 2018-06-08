#include "reaxdetect.h"

int main(int argc, const char** argv) {
	ReaxDetect reaxdetect;
	if (reaxdetect.init(argc, argv))return 1;
	else return reaxdetect.exec();
}