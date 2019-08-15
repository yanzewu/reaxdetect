
#ifndef SERIALIZE_H
#define SERIALIZE_H

#include "external/json.hpp"
#include "argparse.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;




    // Serializer with default arguments
class Serializer {
public:

    Serializer() {

    }

    // get one element
    template<class T>
    T at(const std::string& name)const {
        return _data.at(name).get<T>();
    }

    // get the json object
    const json& get()const {
        return _data;
    }

    // set default
    void set_default(const json& src) {
        _data = src;
    }

    // load file
    void read_file(const std::string& filename) {
        json other;
        std::ifstream infile(filename);
        if (infile.is_open()) {
            infile >> other;
            infile.close();
        }
        else {
            throw std::runtime_error("Cannot open input file");
        }
        merge(other, _data);
    }


    void read_arg(const ArgList& arglist) {
        for (auto d = _data.begin(); d != _data.end(); ++d) {
            if (!d->is_object() && arglist.count(d.key()) > 0) {
                switch (d->type())
                {
                case json::value_t::string:
                    *d = arglist.at(d.key())[0];
                    break;

                case json::value_t::boolean:
                    *d = true;
                    break;

                case json::value_t::number_integer:
                    *d = atoi(arglist.at(d.key())[0].c_str());
                    break;

                case json::value_t::number_unsigned:
                    *d = (unsigned)atoi(arglist.at(d.key())[0].c_str());
                    break;

                case json::value_t::number_float:
                    *d = atof(arglist.at(d.key())[0].c_str());
                    break;

                default:
                    throw std::runtime_error("Cannot read argument");
                    break;
                }
            }
        }
    }

private:

    static void merge(const json& src, json& dst) {
        for (auto src_iter = src.begin(); src_iter != src.end(); src_iter++) {
            if (src_iter->is_object()) {
                merge(*src_iter, dst[src_iter.key()]);
            }
            else {
                dst[src_iter.key()] = *src_iter;
            }
        }
    }

    json _data;
};

#endif
