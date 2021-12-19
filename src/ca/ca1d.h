#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "gateway_key.h"

using namespace std;

// TODO: Reorder template parameters for types... values
// TODO: Supported static type assertions
// TODO: GTest
template<size_t PartitionSize = 3, typename CellType = bool> class CA1D
{
  private:
    vector<CellType> _state;
    GatewayKey<PartitionSize, CellType> _gateway_key;

    void merge_partitions(vector<CellType>& partition, const vector<CellType>& lhs,
                          const vector<CellType>& rhs)
    {
        copy(lhs.begin(), lhs.end(), back_inserter(partition));
        copy(rhs.begin(), rhs.end(), back_inserter(partition));
    }

    /**
     * Retrieve partition with (possibly uneven) radii surrounding the cell
     * @param cell: index of the cell in state
     * @return: standard partition (slice)
     */
    vector<CellType> standard_partition(int cell)
    {
        int rounded_radius = _gateway_key.partition_size() / 2;
        int lhs = cell - rounded_radius;
        int rhs = cell + rounded_radius;

        // Composite output partition from lhs_vec and rhs_vec
        vector<CellType> partition;

        // Adjust for biasing (keeping as a switch for future configurations)
        // FIXME: RHS biasing not working
        if (_gateway_key.partition_size() % 2) {
            switch (_gateway_key.partition_bias()) {
            case PARTITION_BIAS_LHS:
                lhs -= 1;
                break;
            case PARTITION_BIAS_RHS:
                rhs += 1;
                break;
            default:
                break;
            }
        }

        switch (_gateway_key.boundary()) {
        case BOUNDARY_CYCLIC:
            if /* Wrap LHS to end of state */ (lhs < 0) {
                vector<CellType> lhs_vec(_state.begin() + _state.size() + lhs, _state.end());
                vector<CellType> rhs_vec(_state.begin(), _state.begin() + rhs);
                merge_partitions(partition, lhs_vec, rhs_vec);
            } /* Wrap RHS to beginning of state */ else if (rhs > (int) _gateway_key.state_size()) {
                int difference = rhs - _gateway_key.state_size();
                vector<CellType> lhs_vec(_state.begin() + lhs, _state.end());
                vector<CellType> rhs_vec(_state.begin(), _state.begin() + difference);
                merge_partitions(partition, lhs_vec, rhs_vec);
            } /* Otherwise, just slice */ else {
                partition = vector<CellType>(_state.begin() + lhs, _state.begin() + rhs);
            }
            break;
        case BOUNDARY_ZERO:
            if /* Wrap LHS to end of state */ (lhs < 0) {
                vector<CellType> zero_boundary(abs(lhs), 0);
                vector<CellType> in_bounds_state(_state.begin(), _state.begin() + rhs);
                merge_partitions(partition, zero_boundary, in_bounds_state);
            } /* Wrap RHS to beginning of state */ else if (rhs > (int) _gateway_key.state_size()) {
                int difference = rhs - _gateway_key.state_size();
                vector<CellType> in_bounds_state(_state.begin() + lhs, _state.end());
                vector<CellType> zero_boundary(difference, 0);
                merge_partitions(partition, in_bounds_state, zero_boundary);
            } /* Otherwise, just slice */ else {
                partition = vector<CellType>(_state.begin() + lhs, _state.begin() + rhs);
            }
            break;
        default:
            break;
        }

        return partition;
    }

    /**
     * Partition state into a neighborhood or block
     * @param cell: the index of the cell in state
     * @return: partition (slice) adhering to CA type and configuration (i.e., neighborhood, block, etc.)
     */
    vector<CellType> partition(size_t cell)
    {
        switch (_gateway_key.ca_type()) {
        case CA_1D:
            return standard_partition(cell);
        default:
            cerr << "ERROR: CA1D::partition(): unimplemented CA type" << endl;
            exit(-1);
        }
    }

    /**
     * Execute the state transition rule for a single cell
     * @param cell: the index of the cell in state
     * @param rule: the rule number to be used in the interaction
     * @return: new cell state
     */
    CellType cell_interaction(size_t cell, size_t rule)
    {
        vector<CellType> partition = this->partition(cell);

        switch (_gateway_key.interaction()) {
        case INTERACTION_NEIGHBORHOOD_TO_RULE_BIT: {
            // get bit offset of partition in _gateway_key.partition_permutations
            auto permutations = _gateway_key.partition_permutations();
            auto it = find(permutations.begin(), permutations.end(), partition);
            auto bit_offset = distance(permutations.begin(), it);
            return (rule >> bit_offset) & 1;
        }
        case INTERACTION_NEIGHBORHOOD_TO_RULE_BIT_XOR_PREV_CELL:
            break;
        }

        return -1;
    }

  public:
    CA1D(GatewayKey<PartitionSize, CellType> gateway_key)
        : _state(gateway_key.start_state()), _gateway_key(gateway_key)
    {
    }

    /**
     * Evolve the configured CA for @epochs
     * @param epochs: number of timesteps to run evolution for (entire state update)
     */
    void evolve(size_t epochs = 25)
    {
        for (size_t rule = 0; rule < pow(2, _gateway_key.total_permutations()); rule++) {
            cout << "RULE: " << rule << endl;
            vector<vector<bool>> state_history;

            // Reset state from previous runs
            _state = _gateway_key.start_state();
            state_history.push_back(_state);

            for (size_t epoch = 0; epoch < epochs; epoch++) {
                // Capture previous state
                vector<CellType> new_state(_state.begin(), _state.end());

                // Update state
                for (size_t cell = 0; cell < _gateway_key.state_size(); cell++)
                    new_state[cell] = cell_interaction(cell, rule);

                // State transition
                _state = new_state;
                state_history.push_back(_state);
            }

            write_pgm(state_history, rule);
        }
    }

    /**
     * Write state history to image file using state over time
     * @param state_history: 2D state history where rows are states over time
     * @param rule: rule to use as output filename
     */
    void write_pgm(const vector<vector<CellType>>& state_history, int rule)
    {
        // TODO: Create output directory if not exists
        // TODO: Default output directory with configurable option
        ofstream pgm("img/" + to_string(rule) + ".pgm", ios::out | ios::binary);
        // TODO: Assertion that state_history is not empty
        pgm << "P2\n" << state_history[0].size() << " " << state_history.size() << "\n" << 1 << "\n";

        for (auto row : state_history) {
            for (auto pixel : row)
                pgm << pixel << " ";
            pgm << "\n";
        }

        pgm.close();
    }

    GatewayKey<PartitionSize> gateway_key() { return _gateway_key; }

    vector<CellType> state() { return _state; }
};
