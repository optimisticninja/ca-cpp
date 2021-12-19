#include "ca/ca1d.h"

#include "util/states.h"

using namespace std;

/**
 * Terminology:
 *  - partition (synonymous with neighborhood or block)
 *  - radii (or "radius") possibly uneven "midway-point" of partition with biasing for left or right hand side
 * preference
 */

#define DEBUG

int main()
{
    const size_t TOTAL_CELLS = 31;
    const size_t PARTITION_SIZE = 3;
    vector<bool> start_state = wolfram_start_state(TOTAL_CELLS);

    // FIXME: Block CA is very uninteresting, need better block_interaction()
    //     GatewayKey<PARTITION_SIZE, bool> gateway_key(
    //         start_state,
    //         CA_1D_BLOCK,
    //         BOUNDARY_CYCLIC,
    //         PARTITION_BIAS_LHS,
    //         INTERACTION_NEIGHBORHOOD_TO_RULE_BIT_XOR_PREV_CELL
    //     );
    // Standard 1D CA
    GatewayKey<PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_CYCLIC, PARTITION_BIAS_LHS);
    CA1D<PARTITION_SIZE, bool> ca(gateway_key);

    ca.evolve();
}
