#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "../util/special_print.h"
#include "../util/states.h"
#include "gateway_key.h"

using namespace std;

/// Abstract cellular automata
// TODO: Specify dimensionality, lock it down to 1D
template<typename CellType = bool,
         /// Output type of neighborhood (partition) transition
         typename LocalTransitionOutputType = CellType,
         /// Type of state representation
         typename GlobalTransitionOutputType = vector<CellType>,
         /// The number of cells to use when determing state transition for a cell or block
         size_t PartitionSize = 3>
class CA
{
  protected:
    /// State representation
    GlobalTransitionOutputType _state;
    /// Configuration/encoding of the CA
    GatewayKey<PartitionSize, CellType> _gateway_key;

    // Possible FIXME when entering higher dimensions
    /**
     * Merge partitions when reaching the edge of state
     * @param partition: output partition
     * @param lhs: the left hand side of the new partition
     * @param rhs: the right hand side of the new partition
     */
    // TODO: Move into CA1D, this is 1D specific
    void merge_partitions(vector<CellType>& partition, const vector<CellType>& lhs,
                          const vector<CellType>& rhs)
    {
        copy(lhs.begin(), lhs.end(), back_inserter(partition));
        copy(rhs.begin(), rhs.end(), back_inserter(partition));
    }

    virtual vector<CellType> partition(size_t cell) = 0;
    virtual LocalTransitionOutputType local_transition(size_t cell, size_t rule, bool log = true) = 0;
    CA(GatewayKey<PartitionSize, CellType> gateway_key)
        : _state(gateway_key.start_state()), _gateway_key(gateway_key)
    {
    }

  public:
    GatewayKey<PartitionSize> gateway_key() { return _gateway_key; }
    GlobalTransitionOutputType state() { return _state; }
    void state(GlobalTransitionOutputType state) { _state = state; };
    virtual GlobalTransitionOutputType global_transition(size_t epochs = 25, bool write_image = false) = 0;
};

// TODO: Supported static type assertions
// TODO: Create optimized elementary CA for bool only that uses bitsets
// TODO: Condense logging into a single line
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class CA1D : public CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
#ifdef DEBUG
  public:
#endif

    /**
     * Sliding window of operation across entire state. This is the neighborhood size + 1
     * @param cell: index of the target cell in state
     * @return: neighborhood including @cell
     */
    vector<CellType> partition(size_t cell)
    {
        int rounded_radius = this->_gateway_key.partition_size() / 2;
        // FIXME: LHS is a negative index, leftover from cat-playground port
        int lhs = cell - rounded_radius;
        int rhs = cell + rounded_radius;

        // Concatenation of in/out of bounds slices if at edge, otherwise a standard slice
        vector<CellType> partition;

        // Adjust for biasing
        // Possible FIXME: Logic seems off and doesn't need to be switched
        if (!(this->_gateway_key.partition_size() % 2)) {
            switch (this->_gateway_key.partition_bias()) {
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

        // Create partition by merging boundary/slicing
        switch (this->_gateway_key.boundary()) {
        case BOUNDARY_CYCLIC:
            if /* Wrap LHS to end of state */ (lhs < 0) {
                vector<CellType> lhs_vec(this->_state.begin() + this->_state.size() + lhs,
                                         this->_state.end());
                vector<CellType> rhs_vec(this->_state.begin(), this->_state.begin() + rhs);
                this->merge_partitions(partition, lhs_vec, rhs_vec);
            } /* Wrap RHS to beginning of state */ else if ((size_t) rhs > this->_gateway_key.state_size()) {
                int difference = rhs - this->_gateway_key.state_size();
                vector<CellType> lhs_vec(this->_state.begin() + lhs, this->_state.end());
                vector<CellType> rhs_vec(this->_state.begin(), this->_state.begin() + difference);
                this->merge_partitions(partition, lhs_vec, rhs_vec);
            } /* Otherwise, just slice */ else {
                partition = vector<CellType>(this->_state.begin() + lhs, this->_state.begin() + rhs);
            }
            break;
        case BOUNDARY_ZERO:
            if /* Wrap LHS to end of state */ (lhs < 0) {
                vector<CellType> zero_boundary(abs(lhs), 0);
                vector<CellType> in_bounds_state(this->_state.begin(), this->_state.begin() + rhs);
                this->merge_partitions(partition, zero_boundary, in_bounds_state);
            } /* Wrap RHS to beginning of state */ else if ((size_t) rhs > this->_gateway_key.state_size()) {
                int difference = rhs - this->_gateway_key.state_size();
                vector<CellType> in_bounds_state(this->_state.begin() + lhs, this->_state.end());
                vector<CellType> zero_boundary(difference, 0);
                this->merge_partitions(partition, in_bounds_state, zero_boundary);
            } /* Otherwise, just slice */ else {
                partition = vector<CellType>(this->_state.begin() + lhs, this->_state.begin() + rhs);
            }
            break;
        default:
            cerr << "ERROR: CA1D::partition : unsupported boundary type" << endl;
            exit(1);
            break;
        }

        return partition;
    }

    /**
     * Execute the local state transition rule (partition maps to a permutation index used for getting that
     * bit in the rule)
     * @param cell: the index of the cell in state
     * @param rule: the rule number to be used in the interaction
     * @return: new cell state
     */
    LocalTransitionOutputType local_transition(size_t cell, size_t rule, bool log = true)
    {
        auto partition = this->partition(cell);
        cout << "cell:\t\t" << cell << endl;
        cout << "partition:\t";
        print_vector(partition);
        cout << endl;

        switch (this->_gateway_key.interaction()) {
        // Rule Wolfram uses
        case INTERACTION_NEIGHBORHOOD_TO_RULE_BIT: {
            // get bit offset of partition in this->_gateway_key.partition_permutations
            auto permutations = this->_gateway_key.partition_permutations();
            auto it = find(permutations.begin(), permutations.end(), partition);
            auto bit_offset = distance(permutations.begin(), it);
            auto new_state = (rule >> bit_offset) & 1;

            if (log)
                cout << "update:\t\t" << this->_state[cell] << " -> " << new_state << endl;
            return new_state;
        }
        default:
            cerr << "ERROR: CA1D::cell_interaction(): unimplemented interaction" << endl;
            exit(1);
            break;
        }
    }

#ifndef DEBUG
  public:
#endif
    CA1D(GatewayKey<PartitionSize, CellType> gateway_key)
        : CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(gateway_key)
    {
    }

    /**
     * Evolve the configured CA for @epochs (global transition function)
     * @param epochs: number of timesteps to run evolution for (entire state update)
     * @param write_image: write PGM files to observe state over time visually
     * @return: ending state
     */
    // TODO: Should this return state history? Or at all?
    GlobalTransitionOutputType global_transition(size_t epochs = 25, bool write_image = false)
    {
        cout << "start state:\t";
        print_vector(this->gateway_key().start_state());
        cout << endl;

        // Enumerate all rules dependent on the total permutations
        for (size_t rule = 0; rule < pow(2, this->gateway_key().total_permutations()); rule++) {
            cout << "rule:\t" << rule << endl;
            vector<GlobalTransitionOutputType> state_history;

            // Reset state to initial state
            this->_state = this->gateway_key().start_state();
            state_history.push_back(this->_state);

            for (size_t epoch = 0; epoch < epochs; epoch++) {
                // Capture previous state
                GlobalTransitionOutputType new_state(this->_state.begin(), this->_state.end());

                for (size_t cell = 0; cell < new_state.size(); cell++)
                    new_state[cell] = local_transition(cell, rule);

                // State transition
                this->state(new_state);
                state_history.push_back(this->state());
            }

            if (write_image)
                write_pgm(state_history, rule);
        }

        return this->_state;
    }
};
