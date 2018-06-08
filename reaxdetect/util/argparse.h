#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <string>
#include <vector>
#include <stdexcept>
#include <ostream>
#include <map>

class ArgParseError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

typedef std::map<std::string, std::vector<std::string> > ArgList;

class ArgParser {
public:

    struct Option {
        std::string name;
        std::string short_name;
        std::string long_name;
        unsigned nargs;
    };

    ArgParser() : _show_help(true) {
        options.push_back({ "help", "-h", "--help", 0});
    }

    ArgParser(bool show_help) : _show_help(show_help) {
        if (_show_help) {
            options.push_back({ "help", "-h", "--help", 0});
        }
    }

    // require_arg: 0/1/+
    void add_argument(const std::string& short_name, const std::string& long_name, char nargs = 1) {
        if (!_is_option(long_name) || !_is_option(short_name) || !_is_long_option(long_name)) {
            throw std::runtime_error("Not a valid option");
        }
        options.push_back({ long_name.substr(2), short_name, long_name, nargs == '+' ? (unsigned)-1 : nargs});
    }

    void add_argument(const std::string& name, char nargs = 1) {
        if (_is_option(name)) {
            if (!_is_long_option(name)) {
                throw std::runtime_error("Not a long option or positional option");
            }
            options.push_back({ name.substr(2), "", name, nargs == '+' ? (unsigned)-1 : nargs});
        }
        else {
            if (nargs == 0)nargs = 1;
            positional_options.push_back({ name, "", "", nargs == '+' ? (unsigned)-1 : nargs});
        }
    }

    void set_prog_name(const std::string& name) {
        prog_name = name;
    }

    // first argument will be ignored
    ArgList parse(int argc, const char** argv) {
        std::vector<std::string> strlist;

        for (int i = 1; i < argc; i++) {
            strlist.push_back(argv[i]);
        }
        return parse(strlist);
    }

    ArgList parse(const std::vector<std::string>& strlist) {
        ArgList result;
        auto parg = positional_options.begin();

        for (auto s = strlist.begin(); s != strlist.end();) {
            if (!_is_option(*s)) {
                if (parg != positional_options.end()){
                    for (unsigned i = 0; i < parg->nargs; ++i) {
                        if (s == strlist.end() || _is_option(*s)) {
                            throw ArgParseError("Not enough args: " + parg->name);
                        }
                        else {
                            result[parg->name].push_back(*s);
                            ++s;
                        }
                    }
                    ++parg;
                }
                else {
                    throw ArgParseError("Unrecognized positional arg: " + *s);
                }
            }
            else {
                const auto& opt = loc_option(*s);
                ++s;

                result[opt.name] = std::vector<std::string>();
                for (unsigned i = 0; i < opt.nargs; ++i) {
                    if (s == strlist.end() || _is_option(*s)) {
                        throw ArgParseError("Not enough args: " + opt.name);
                    }
                    else {
                        result[opt.name].push_back(*s);
                        ++s;
                    }
                }
            }
        }
        return result;
    }

    void display_help(std::ostream& os) {
        
        auto to_upper = [](const std::string& s){
            std::string ret;
            for (const auto c : s)ret.push_back(toupper(c));
            return ret;
        };

        os << "Usage: ";
        os << prog_name;
        for (const auto& opt : options) {
            os << " [";
            if (!opt.short_name.empty()) {
                os << opt.short_name;
            }
            else {
                os << opt.long_name;
            }
            if (opt.nargs >= 1) {
                os << ' ' << to_upper(opt.name);
            }
            if (opt.nargs > 1) {
                os << " ...";
            }
            os << ']';
        }
        for (const auto& popt : positional_options) {
            os << ' ' << to_upper(popt.name);
            if (popt.nargs >= 1) {
                os << " ...";
            }
        }
        os << std::endl;
    }

private:

    static bool _is_option(const std::string& str) {
        return str.length() > 1 && str[0] == '-';
    }

    static bool _is_long_option(const std::string& str) {
        return str.length() > 2 && str[0] == '-' && str[1] == '-';
    }

    // if not found ==> throw error
    const Option& loc_option(const std::string& name) {
        if (_is_long_option(name)) {
            for (const auto& opt : options) {
                if (opt.long_name == name) {
                    return opt;
                }
            }
        }
        else {
            for (const auto& opt : options) {
                if (opt.short_name == name) {
                    return opt;
                }
            }
        }
        throw ArgParseError("Invalid argument: " + name);
    }

    bool _show_help;
    std::string prog_name;
    std::vector<Option> options;
    std::vector<Option> positional_options;
};

#endif