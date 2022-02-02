#include "states.h"

#include <random>

using namespace std;

static random_device rd;
static mt19937 mt(rd());

vector<bool> wolfram_start_state(size_t total_cells)
{
    vector<bool> start_state;
    for (size_t i = 0; i < total_cells; i++)
        start_state.push_back(true);
    start_state[total_cells / 2] = false;
    return start_state;
}

vector<vector<bool>> glider(size_t x, size_t y)
{
    vector<vector<bool>> start_state(y);

    for (size_t i = 0; i < y; i++)
        for (size_t j = 0; j < x; j++)
            start_state[i].push_back(0);

    start_state[1][2] = true;
    start_state[2][3] = true;
    start_state[3][1] = true;
    start_state[3][2] = true;
    start_state[3][3] = true;

    return start_state;
}

vector<bool> random_1d_start_state(size_t size, double probability)
{
    vector<bool> start_state(size);
    bernoulli_distribution dist(probability);
    for (size_t i = 0; i < size; i++)
        start_state[i] = dist(mt);
    return start_state;
}

vector<vector<bool>> random_2d_start_state(size_t x, size_t y, double probability)
{
    vector<vector<bool>> start_state(y);
    bernoulli_distribution dist(probability);

    for (size_t i = 0; i < y; i++) {
        vector<bool> row_state(x);
        for (size_t j = 0; j < x; j++)
            row_state[j] = dist(mt);
        start_state[i] = row_state;
    }

    return start_state;
}
