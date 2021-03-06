// #include <bits/stdc++.h>
#include <algorithm>

#include "ca/ca1d.h"
#include "ca/ca2d.h"

#include "util/states.h"

using namespace std;

#define DEBUG

const size_t STATE_SIZE = 480;
const size_t EPOCHS = 360;
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
    vector<vector<bool>> start_state = random_2d_start_state(STATE_SIZE, STATE_SIZE);
    GatewayKey2D<vector<vector<bool>>, 8, bool> gateway_key(start_state, CA_2D);
    Alias2D::Life ca(gateway_key);
    ca.evolve(10000, WRITE_IMAGES);
}

void second_order_ca()
{
    vector<bool> start_state = random_1d_start_state(STATE_SIZE);
    vector<bool> prev_state = random_1d_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>> gateway_key(start_state, prev_state, CA_1D, BOUNDARY_CYCLIC,
                                           PARTITION_BIAS_LHS,
                                           INTERACTION_SECOND_ORDER_NEIGHBORHOOD_TO_RULE_BIT);
    Alias1D::SecondOrderCA ca(gateway_key);
    //     ca.evolve_all(EPOCHS, WRITE_IMAGES);
    vector<vector<bool>> forward = ca.evolve_rule(90, EPOCHS);
    vector<vector<bool>> backward =
        ca.evolve_rule_reverse(90, EPOCHS, forward[forward.size() - 1], forward[forward.size() - 2]);
    vector<vector<bool>> backward_reversed(backward.rbegin(), backward.rend());
    write_pgm(forward, 90, "reversed");
    write_pgm(backward_reversed, 90, "reversed", "R");
}

int main()
{
    //     elementary_ca();
    second_order_ca();
    //     game_of_life();
}
