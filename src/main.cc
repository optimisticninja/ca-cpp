#include "ca/ca1d.h"

using namespace std;

/**
 * Terminology:
 *  - partition (synonymous with neighborhood or block)
 *  - radii (or "radius") possibly uneven "midway-point" of partition with biasing for left or right hand side
 * preference
 */

vector<bool> wolfram_start_state(size_t total_cells)
{
    vector<bool> start_state;
    for (size_t i = 0; i < total_cells; i++)
        start_state.push_back(true);
    start_state[total_cells / 2 + 1] = false;
    return start_state;
}

int main()
{
    const size_t TOTAL_CELLS = 31;
    const size_t PARTITION_SIZE = 3;
    vector<bool> start_state = wolfram_start_state(TOTAL_CELLS);

    // Configure the CA
    GatewayKey<PARTITION_SIZE, bool> gateway_key(start_state, BOUNDARY_ZERO, PARTITION_BIAS_RHS);
    CA1D<PARTITION_SIZE, bool> ca(gateway_key);

    ca.evolve();
}
