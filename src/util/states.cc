#include "states.h"

#include <random>

vector<bool> wolfram_start_state(size_t total_cells)
{
    vector<bool> start_state;
    for (size_t i = 0; i < total_cells; i++)
        start_state.push_back(true);
    start_state[total_cells / 2] = false;
    return start_state;
}

vector<vector<bool>> random_2d_start_state(size_t x, size_t y)
{
    vector<vector<bool>> start_state(y);

    for (size_t i = 0; i < y; i++)
        for (size_t j = 0; j < x; j++)
            start_state[i].push_back(rand() % 2);

    return start_state;
}
