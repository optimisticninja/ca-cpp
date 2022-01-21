#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

using namespace std;

vector<bool> wolfram_start_state(size_t total_cells);

vector<vector<bool>> random_2d_start_state(size_t x, size_t y);
vector<vector<bool>> glider(size_t x, size_t y);

/**
 * Write state history to image file using state over time
 * @param state_history: 2D state history where rows are states over time
 * @param rule: rule to use as output filename
 */
// TODO: Use compressed image format (preferably PNG)
template<typename GlobalTransitionOutputType>
void write_pgm(const vector<GlobalTransitionOutputType>& state_history, int rule,
               const char* directory = "img", const char* suffix = "")
{
    if (state_history.empty()) {
        cerr << "ERROR: CA1D::write_pgm(): empty state history" << endl;
        exit(1);
    }

    // Create directory if doesn't exist
    // TODO: Configurable output directory
    struct stat st;
    if (stat(directory, &st) == -1) {
        // TODO: Correct permissions
        mkdir(directory, 0700);
    }

    ofstream pgm(string(directory) + "/" + to_string(rule) + string(suffix) + ".pgm", ios::out | ios::binary);
    pgm << "P2\n" << state_history[0].size() << " " << state_history.size() << "\n" << 1 << "\n";

    for (auto row : state_history) {
        for (auto pixel : row)
            pgm << pixel << " ";
        pgm << "\n";
    }

    pgm.close();
};

/**
 * Write 2D state timestep to image file
 * @param state_history: 2D state history where rows are states over time
 * @param rule: rule to use as output filename
 */
// TODO: Use compressed image format (preferably PNG)
template<typename GlobalTransitionOutputType>
void write_pgm_2d_state(const GlobalTransitionOutputType& timestep, int epoch, const char* directory = "img")
{
    if (timestep.empty()) {
        cerr << "ERROR: CA1D::write_pgm(): empty state history" << endl;
        exit(1);
    }

    // Create directory if doesn't exist
    // TODO: Configurable output directory
    struct stat st;
    if (stat(directory, &st) == -1) {
        // TODO: Correct permissions
        mkdir(directory, 0700);
    }

    ostringstream ss;
    ss << setw(6) << setfill('0') << epoch;
    ofstream pgm(string(directory) + "/" + ss.str() + ".pgm", ios::out | ios::binary);
    pgm << "P2\n" << timestep[0].size() << " " << timestep.size() << "\n" << 1 << "\n";

    for (auto row : timestep) {
        for (auto pixel : row)
            pgm << pixel << " ";
        pgm << "\n";
    }

    pgm.close();
};
