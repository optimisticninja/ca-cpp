#include "ca/ca1d.h"

#include "util/states.h"

using namespace std;

#define DEBUG

void elementary_ca()
{
    const size_t STATE_SIZE = 31;
    const size_t EPOCHS = 25;
    const bool WRITE_IMAGES = true;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey<> gateway_key(start_state, CA_1D, BOUNDARY_CYCLIC);
    Alias1D::ElementaryCA ca(gateway_key);
    ca.evolve_all(EPOCHS, WRITE_IMAGES);
}

void game_of_life()
{
    // TODO: Pass dimensionality as vector of shape (i.e. { 2, 2 }) to gateway key as means of scaling to
    // higher dimensions
}

int main() { elementary_ca(); }
