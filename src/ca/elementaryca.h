#pragma once

#include <bitset>
#include <cassert>
#include <vector>

using namespace std;

typedef enum { BOUNDARY_PERIODIC, BOUNDARY_ZERO } boundary_t;

/** Optimized elementary CA using bitsets */
// TODO: Align to template interfaces, this was in use for a novel encryption algorithm and should be adapted
template<size_t StateSize> class ElementaryCA
{
  private:
    bitset<StateSize> _state;
    vector<bitset<3>> permutations;
    boundary_t boundary;

    /**
     * Recursive function to generate all permutations of specified neighborhood cells. Used
     * for finding the index of a partition within the exhaustive set and using it as a
     * bit-shift within qualifying rules.
     * @param permutations: vector to store generated permutations in
     * @param neighborhood: starting neighborhood to permute
     * @param cell: cell to operate on
     */
    void generate_partition_states(vector<bitset<3>>& permutations, bitset<3>& neighborhood, size_t cell)
    {
        if (cell == 3) {
            permutations.push_back(neighborhood);
            return;
        }

        neighborhood[cell] = 0;
        generate_partition_states(permutations, neighborhood, cell + 1);
        neighborhood[cell] = 1;
        generate_partition_states(permutations, neighborhood, cell + 1);
    }

  public:
    ElementaryCA(const bitset<StateSize>& state, boundary_t boundary = BOUNDARY_PERIODIC)
        : _state(state), boundary(boundary)
    {
        bitset<3> neighborhood(3);
        generate_partition_states(permutations, neighborhood, 0);
    }

    /**
     * Run CA for 1 epoch
     * @param rule CA transition rule to execute
     */
    bitset<StateSize> step(unsigned rule)
    {
        bitset<StateSize> new_state(_state);

        if (boundary == BOUNDARY_PERIODIC) {
            for (unsigned cell = 0; cell < _state.size(); cell++) {
                bitset<3> neighborhood;

                // Create a neighborhood with periodic boundary
                if (cell == 0) {
                    neighborhood[0] = _state[_state.size() - 1];
                    neighborhood[1] = _state[cell];
                    neighborhood[2] = _state[cell + 1];
                } else if (cell == (_state.size() - 1)) {
                    neighborhood[0] = _state[cell - 1];
                    neighborhood[1] = _state[cell];
                    neighborhood[2] = _state[0];
                } else {
                    neighborhood[0] = _state[cell - 1];
                    neighborhood[1] = _state[cell];
                    neighborhood[2] = _state[cell + 1];
                }

                // get bit offset of partition in this->_gateway_key.partition_permutations
                unsigned i = 0;
                for (; i < permutations.size(); i++)
                    if (permutations[i] == neighborhood)
                        break;
                new_state[cell] = (rule >> i) & 1;
            }
        } else if (boundary == BOUNDARY_ZERO) {
            for (unsigned cell = 0; cell < _state.size(); cell++) {
                bitset<3> neighborhood;

                // Create a neighborhood with null boundary
                if (cell == 0) {
                    neighborhood[0] = 0;
                    neighborhood[1] = _state[cell];
                    neighborhood[2] = _state[cell + 1];
                } else if (cell == (_state.size() - 1)) {
                    neighborhood[0] = _state[cell - 1];
                    neighborhood[1] = _state[cell];
                    neighborhood[2] = 0;
                } else {
                    neighborhood[0] = _state[cell - 1];
                    neighborhood[1] = _state[cell];
                    neighborhood[2] = _state[cell + 1];
                }

                // get bit offset of neighborhood in permutations
                unsigned i = 0;
                for (; i < permutations.size(); i++)
                    if (permutations[i] == neighborhood)
                        break;
                new_state[cell] = (rule >> i) & 1;
            }
        } else {
            assert(true && "Unsupported boundary type");
        }

        this->_state = new_state;
        return new_state;
    }

    bitset<StateSize> state() { return _state; }
};
