#include "ca/ca1d.h"
#include "ca/ca2d.h"

#include "util/states.h"

using namespace std;

#define DEBUG

const size_t STATE_SIZE = 31;
const size_t EPOCHS = 25;
const bool WRITE_IMAGES = true;

void elementary_ca()
{
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>> gateway_key(start_state, CA_1D, BOUNDARY_CYCLIC);
    Alias1D::ElementaryCA ca(gateway_key);
    ca.evolve_all(EPOCHS, WRITE_IMAGES);
}

void game_of_life()
{
    vector<vector<bool>> start_state = random_2d_start_state(31, 31);
    GatewayKey2D<vector<vector<bool>>, 8, bool> gateway_key(start_state, CA_2D, BOUNDARY_CYCLIC);
    Alias2D::Life ca(gateway_key);
    // TODO: Pass dimensionality as vector of shape (i.e. { 2, 2 }) to gateway key as means of scaling to
    // higher dimensions
}

int main() { elementary_ca(); }
