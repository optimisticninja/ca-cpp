
#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "../util/special_print.h"
#include "gateway_key.h"

using namespace std;

// TODO: Reorder template parameters for types... values
// TODO: Supported static type assertions
template<size_t PartitionSize = 3, typename CellType = bool> class CA1D
{
  private:
    vector<CellType> _state;
    GatewayKey<PartitionSize, CellType> _gateway_key;

#ifdef DEBUG
  public:
#endif
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
        // FIXME: LHS is a negative index, leftover from cat-playground port
        int lhs = cell - rounded_radius;
        int rhs = cell + rounded_radius;

        // Concatenation of bounded slices from lhs/rhs if at edge, otherwise a standard slice from array
        vector<CellType> partition;

        // Adjust for biasing
        if (!(_gateway_key.partition_size() % 2)) {
            switch (_gateway_key.partition_bias()) {
            case PARTITION_BIAS_RHS:
                // Shift cell partition over 1
                rhs += 1;
                lhs += 1;
                break;
            default:
                break;
            }
        } else {
            rhs += 1;
        }

        switch (_gateway_key.boundary()) {
        case BOUNDARY_CYCLIC:
            if /* Wrap LHS to end of state */ (lhs < 0) {
                vector<CellType> lhs_vec(_state.begin() + _state.size() + lhs, _state.end());
                vector<CellType> rhs_vec(_state.begin(), _state.begin() + rhs);
                merge_partitions(partition, lhs_vec, rhs_vec);
            } /* Wrap RHS to beginning of state */ else if ((size_t) rhs > _gateway_key.state_size()) {
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
            } /* Wrap RHS to beginning of state */ else if ((size_t) rhs > _gateway_key.state_size()) {
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

    vector<CellType> block_partition(size_t cell)
    {
        // Concatenation of bounded slices from lhs/rhs if at edge, otherwise a standard slice from array
        vector<CellType> partition;

        // TODO: Adjust for biasing
        if (!(_gateway_key.partition_size() % 2)) {
            switch (_gateway_key.partition_bias()) {
            case PARTITION_BIAS_RHS:
                break;
            case PARTITION_BIAS_CENTER:
                break;
            default:
                break;
            }
        }

        size_t difference = _gateway_key.state_size() - cell;
        int bounded_cells = difference - _gateway_key.partition_size();

        switch (_gateway_key.boundary()) {
        case BOUNDARY_CYCLIC:
            if (bounded_cells < 0) {
                vector<CellType> boundary(_state.begin(), _state.begin() + abs(bounded_cells));
                vector<CellType> in_bounds(_state.begin() + cell, _state.end());
                merge_partitions(partition, in_bounds, boundary);
            } else {
                partition = vector<CellType>(_state.begin() + cell,
                                             _state.begin() + cell + _gateway_key.partition_size());
            }
            break;
        case BOUNDARY_ZERO:
            if (bounded_cells < 0) {
                vector<CellType> boundary(abs(bounded_cells), 0);
                vector<CellType> in_bounds(_state.begin() + cell, _state.end());
                merge_partitions(partition, in_bounds, boundary);
            } else {
                partition = vector<CellType>(_state.begin() + cell,
                                             _state.begin() + cell + _gateway_key.partition_size());
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
        case CA_1D_BLOCK:
            return block_partition(cell);
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
    CellType cell_interaction(size_t cell, size_t rule, bool log = true)
    {
        auto partition = this->partition(cell);
        cout << "cell:\t\t" << cell << endl;
        cout << "partition:\t";
        print_vector(partition);
        cout << endl;

        switch (_gateway_key.interaction()) {
        // Rule Wolfram uses
        case INTERACTION_NEIGHBORHOOD_TO_RULE_BIT: {
            // get bit offset of partition in _gateway_key.partition_permutations
            auto permutations = _gateway_key.partition_permutations();
            auto it = find(permutations.begin(), permutations.end(), partition);
            auto bit_offset = distance(permutations.begin(), it);
            auto new_state = (rule >> bit_offset) & 1;

            if (log)
                cout << "update:\t\t" << _state[cell] << " -> " << new_state << endl;
            return new_state;
        }
        default:
            cerr << "ERROR: CA1D::cell_interaction(): unimplemented interaction" << endl;
            exit(1);
            break;
        }
    }

    /**
     * Execute the state transition rule for a block
     * @param cell: the index of the cell in state
     * @param rule: the rule number to be used in the interaction
     * @return: new cell state
     */
    vector<CellType> block_interaction(size_t cell, size_t rule, bool log = true)
    {
        auto partition = this->partition(cell);
        cout << "partition:\t";
        print_vector(partition);
        cout << endl;

        // TODO: Trim excess from block added during translation
        switch (_gateway_key.interaction()) {
        case INTERACTION_NEIGHBORHOOD_TO_RULE_BIT_XOR_PREV_CELL: {
            // TODO: Restrict to boolean CAs
            // get bit offset of partition in _gateway_key.partition_permutations
            auto permutations = _gateway_key.partition_permutations();
            auto it = find(permutations.begin(), permutations.end(), partition);
            auto bit_offset = distance(permutations.begin(), it);

            // For logging only
            vector<CellType> old_state;
            vector<CellType> new_state;

            for (size_t i = 0; i < _gateway_key.partition_size(); i++) {
                CellType prev_cell;
                old_state.push_back(_state.at(cell));

                // NOTE: Implicitly cyclic and lhs biased
                prev_cell = (cell == 0) ? _state.at(_state.size() - 1) : _state.at(cell - 1);

                // FIXME: Is this second-order in a block since we are XORing?
                CellType new_cell_state = (rule >> bit_offset) ^ _state[cell - 1] & 1;
                new_state.push_back(new_cell_state);
            }

            if (log) {
                cout << "update:\t\t";
                print_vector(old_state);
                cout << " -> ";
                print_vector(new_state);
                cout << endl;
            }

            return new_state;
        }
        default:
            cerr << "ERROR: CA1D::block_interaction() : unsupported interaction" << endl;
            exit(-1);
        }
    }

#ifndef DEBUG
  public:
#endif
    CA1D(GatewayKey<PartitionSize, CellType> gateway_key)
        : _state(gateway_key.start_state()), _gateway_key(gateway_key)
    {
    }

    /**
     * Evolve the configured CA for @epochs
     * @param epochs: number of timesteps to run evolution for (entire state update)
     */
    void evolve(size_t epochs = 25, bool write_image = false)
    {
        switch (_gateway_key.ca_type()) {
        case CA_1D:
            for (size_t rule = 0; rule < pow(2, _gateway_key.total_permutations()); rule++) {
                cout << "rule:\t" << rule << endl;
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

                if (write_image)
                    write_pgm(state_history, rule);
            }
            break;
        case CA_1D_BLOCK:
            for (size_t rule = 0; rule < pow(2, _gateway_key.total_permutations()); rule++) {
                cout << "rule:\t" << rule << endl;
                vector<vector<bool>> state_history;

                // Reset state from previous runs
                _state = _gateway_key.start_state();
                state_history.push_back(_state);

                for (size_t epoch = 0; epoch < epochs; epoch++) {
                    for (size_t block_start = 0; block_start < _state.size();
                         block_start += _gateway_key.partition_size()) {
                        vector<CellType> new_state = block_interaction(block_start, rule);

                        // TEST: Splice new_state into _state
                        for (size_t i = block_start; i < _gateway_key.partition_size(); i++)
                            _state[i] = new_state[i];
                    }

                    state_history.push_back(_state);
                }

                if (write_image)
                    write_pgm(state_history, rule);
            }
            break;
        default:
            break;
        }
    }

    /**
     * Write state history to image file using state over time
     * @param state_history: 2D state history where rows are states over time
     * @param rule: rule to use as output filename
     */
    // TODO: Use compressed image format (preferably PNG)
    void write_pgm(const vector<vector<CellType>>& state_history, int rule, const char* directory = "img")
    {
        if (state_history.empty()) {
            cerr << "ERROR: CA1D::write_pgm(): empty state history" << endl;
            exit(1);
        }

        // Create directory if doesn't exist
        // TODO: Configurable output directory
        struct stat st;
        if (stat(directory, &st) == -1) {
            // TODO: Correct permissions
            mkdir(directory, 0700);
        }

        ofstream pgm(string(directory) + "/" + to_string(rule) + ".pgm", ios::out | ios::binary);
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
