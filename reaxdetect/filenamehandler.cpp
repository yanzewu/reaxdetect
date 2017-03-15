#include "filenamehandler.h"
#include "stringconvert.h"

#ifdef __unix__
#include <unistd.h>
#else
#include <io.h>
#define access _access
#define F_OK 0
#endif // __unix__

namespace path {
	bool exists(const string& path) {
		return access(path.c_str(), F_OK) == 0;
	}

	string get_path(const string& path) {
		if (!exists(path)) {
			return path;
		}

		long long filecnt = 2;
		auto path_split = split(path);
		string result;
		while(exists(result = path_split[0] + path_split[1] + "_" + to_string(filecnt) + "." + path_split[2])){
			filecnt++;
		}
		return result;
	}

	vector<string> split(const string& pathIn) {

		vector<string> out(3);
		size_t line = pathIn.find_last_of('/');
		size_t point = pathIn.find_last_of('.');
		if (line == string::npos)line = pathIn.find_last_of('\\');
		if (line == string::npos)out[0] = "";
		else out[0] = pathIn.substr(0, line + 1);
		out[1] = pathIn.substr(line + 1, point - line - 1);
		out[2] = pathIn.substr(point + 1, pathIn.length() - point);
		return out;
	}
}


