#include "states.h"

vector<bool> wolfram_start_state(size_t total_cells)
{
    vector<bool> start_state;
    for (size_t i = 0; i < total_cells; i++)
        start_state.push_back(true);
    start_state[total_cells / 2] = false;
    return start_state;
}
