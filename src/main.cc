#include "ca/ca1d.h"
#include "ca/ca2d.h"

#include "util/states.h"

using namespace std;

#define DEBUG

const size_t STATE_SIZE = 31;
const size_t EPOCHS = 1000;
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
    vector<vector<bool>> start_state = random_2d_start_state(64, 64);
    GatewayKey2D<vector<vector<bool>>, 8, bool> gateway_key(start_state, CA_2D);
    Alias2D::Life ca(gateway_key);
    ca.evolve(EPOCHS, WRITE_IMAGES);
}

int main()
{ /*elementary_ca();*/
    srand(time(NULL));
    game_of_life();
}
