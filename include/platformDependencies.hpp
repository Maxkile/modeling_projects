#ifndef PLATFORM_HPP
#define PLATFORM_HPP

#pragma once

#ifdef __linux__

#include <sys/stat.h>

using OPENMP_INDEX_TYPE = size_t;

#elif _WIN64

#include <direct.h>
#define mkdir(filename, mode) (_mkdir(filename))
using OPENMP_INDEX_TYPE = long long;

#elif _WIN32

#include <direct.h>
#define mkdir(filename, mode) (_mkdir(filename))
using OPENMP_INDEX_TYPE = long long;

#elif __APPLE__

#include <sys/stat.h>

using OPENMP_INDEX_TYPE = size_t;

#endif

#endif
