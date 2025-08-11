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

void create_output_directory() {
    // Create the "output" directory if it doesn't exist
    const char* output_directory = "output";
    struct stat info;
    if (stat(output_directory, &info) != 0) {
        int status = mkdir(output_directory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (status == 0) {
            cout << "Directory '" << output_directory << "' created successfully." << endl;
        } else {
            cerr << "Failed to create the directory '" << output_directory << "'" << endl;
            exit(-1);
        }
    } else {
        cout << "Directory '" << output_directory << "' already exists." << endl;
    }
}

void printResult(string label, double y, double pt, double result, double error) {
    cout.width(3);
    cout << label;
    cout.width(5);
    cout << y;
    cout.width(5);
    cout << pt;
    cout.width(10);
    cout << result;
    cout.width(13);
    cout << error << endl;
}

