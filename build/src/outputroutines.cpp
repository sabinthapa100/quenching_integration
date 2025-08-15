/*
 
 outputroutines.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <fcntl.h>
#include <limits>

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#include "outputroutines.h"

using namespace std;

void print_line() {
    for (int i=0;i<100;i++) cout << "-"; cout << endl;
    return;
}

static bool ensure_dir(const char* path) {
    struct stat info;
    if (stat(path, &info) == 0) return S_ISDIR(info.st_mode);
    if (mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) return true;
    cerr << "Failed to create directory '" << path << "': " << strerror(errno) << "\n";
    return false;
}
void create_output_directory() {
    // Create the "output" directory if it doesn't exist
    const char* baseDir = "output";
    struct stat info;
    if (stat(baseDir, &info) != 0) {
        if (mkdir(baseDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) {
            cout << "Directory '" << baseDir << "' created successfully." << endl;
        } else {
            cerr << "Failed to create the directory '" << baseDir << "': "
                 << strerror(errno) << endl;
            exit(-1);
        }
    } else if (!S_ISDIR(info.st_mode)) {
        cerr << "'" << baseDir << "' exists but is not a directory." << endl;
        exit(-1);
    } else {
        // cout << "Directory '" << baseDir << "' already exists." << endl;
    }
}

void create_output_directory(const std::string& particleName) {
    // Base output folder
    const char* baseDir = "output";
    struct stat info;

    // Create base output folder if it doesn't exist
    if (stat(baseDir, &info) != 0) {
        if (mkdir(baseDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) {
            cout << "Directory '" << baseDir << "' created successfully." << endl;
        } else {
            cerr << "Failed to create directory '" << baseDir << "'" << endl;
            exit(-1);
        }
    }

    // Create particle-specific subfolder
    std::string particleDir = std::string(baseDir) + "/" + particleName;
    if (stat(particleDir.c_str(), &info) != 0) {
        if (mkdir(particleDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) {
            cout << "Directory '" << particleDir << "' created successfully." << endl;
        } else {
            cerr << "Failed to create directory '" << particleDir << "'" << endl;
            exit(-1);
        }
    }
}


void printResult(string label, double y, double pt, double result, double error) {
    if (!std::isfinite(result)) result = 0.0;
    if (!std::isfinite(error))  error  = 0.0;
    cout.width(3);
    cout << label << endl;
    cout.width(5);
    cout << y;
    cout.width(5);
    cout << pt;
    cout.width(10);
    cout << result;
    cout.width(13);
    cout << error << endl;
}

