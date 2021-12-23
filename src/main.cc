#include "ca/ca1d.h"

#include "util/states.h"

using namespace std;

#define DEBUG

int main()
{
    const size_t TOTAL_CELLS = 31;
    const size_t EPOCHS = 25;
    const bool WRITE_IMAGES = true;

    // Configure CA
    vector<bool> start_state = wolfram_start_state(TOTAL_CELLS);
    GatewayKey<> gateway_key(start_state, CA_1D, BOUNDARY_CYCLIC);

    // Create CA
    //  - template parameters are cell type, local transition output type, global transition output type and
    //  partition size
    CA1D<> ca(gateway_key);

    ca.global_transition(EPOCHS, WRITE_IMAGES);
}
