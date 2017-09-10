#pragma once

#ifndef ERRORS_H
#define ERRORS_H

#include <stdexcept>
#include <string>

class ReaxDetectError : public std::runtime_error {
public:

    ReaxDetectError() : std::runtime_error("Unknown Error") {

    }
    ReaxDetectError(char* const str) : std::runtime_error(str) {

    }
    ReaxDetectError(const std::string& str) : std::runtime_error(str) {

    }
};

class IOError : public ReaxDetectError {
public:

    IOError(const std::string& filename) : ReaxDetectError("Cannot open " + filename) {

    }
};

#endif