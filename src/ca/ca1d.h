#pragma once

#include <algorithm>
#include <cmath>
// #include <concepts>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "../util/special_print.h"
#include "../util/states.h"
#include "gateway_key.h"

// TODO: Refactor using the word "partition", apparently there are PCAs, so this term is not a good generic
// TODO: Genericize partition per dimension (use attributes to create blocks or neighborhoods)
// TODO: Use 'concept' from c++20 to restrict to supported types
// TODO: Evolve rule, all doesn't scale to specific 2D implementations like life. Remove from base class
//     concept Derived = std::is_same<U, Class1>::value || std::is_same<U, Class2>::value;
using namespace std;

/**
 * @brief abstract cellular automata
 * @tparam LocalTransitionOutputType the output type of the local transition function (cellular primitive or
 * collection for block output)
 * @tparam GlobalTransitionOutputType the output type of the global transition function (type of state
 * representation)
 * @tparam PartitionSize the size of the sliding window used to iterate over state for a transition (including
 * target cell)
 */
// TODO: Derive GlobalTransitionOutputType after GateWay key accepts shape of CA as a vector to be
// multidimensional
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class CA
{
  protected:
    /// State representation
    GlobalTransitionOutputType _state;
    /// Configuration/encoding of the CA
    GatewayKey<PartitionSize, CellType> _gateway_key;

    CA(GatewayKey<PartitionSize, CellType> gateway_key)
        : _state(gateway_key.start_state()), _gateway_key(gateway_key)
    {
    }

  public:
    GatewayKey<PartitionSize> gateway_key() { return _gateway_key; }
    GlobalTransitionOutputType state() { return _state; }
    void state(GlobalTransitionOutputType state) { _state = state; };
};

template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class CA1D : public CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
  protected:
    /**
     * @brief merge partitions when reaching the edge of state
     * @param partition output partition
     * @param lhs the left hand side of the new partition
     * @param rhs the right hand side of the new partition
     */
    void merge_partitions(vector<CellType>& partition, const vector<CellType>& lhs,
                          const vector<CellType>& rhs)
    {
        // TODO: In-place update of pre-allocated vector is probably faster
        copy(lhs.begin(), lhs.end(), back_inserter(partition));
        copy(rhs.begin(), rhs.end(), back_inserter(partition));
    }

    CA1D(GatewayKey<PartitionSize, CellType> gateway_key)
        : CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(gateway_key)
    {
    }
};

// TODO: Create optimized elementary CA for bool only that uses bitsets
// TODO: Condense logging into a single line
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class IrreversibleCA1D
    : public CA1D<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
#ifdef DEBUG
  public:
#endif

    /**
     * Sliding window of operation across entire state. This is the neighborhood size + 1
     * @param cell index of the target cell in state
     * @return neighborhood including @cell
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
     * @param cell the index of the cell in state
     * @param rule the rule number to be used in the interaction
     * @param log log single-cell state changes
     * @return: new cell state
     */
    LocalTransitionOutputType local_transition(size_t cell, size_t rule, bool log = false)
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
            cerr << "ERROR: CA1D::local_transition(): unimplemented transition" << endl;
            exit(1);
            break;
        }
    }

    /**
     * Evolve the full state from one timestep to another
     * @param rule rule number to use in local transition
     * @return ending state
     */
    GlobalTransitionOutputType global_transition(size_t rule)
    {
        cout << "start state:\t";
        print_vector(this->gateway_key().start_state());
        cout << endl;

        GlobalTransitionOutputType new_state(this->_state.begin(), this->_state.end());

        for (size_t cell = 0; cell < new_state.size(); cell++)
            new_state[cell] = local_transition(cell, rule);

        // State transition
        this->state(new_state);
        cout << "state: ";
        print_vector(this->_state);
        cout << endl;
        return this->_state;
    }

#ifndef DEBUG
  public:
#endif
    IrreversibleCA1D(GatewayKey<PartitionSize, CellType> gateway_key)
        : CA1D<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(gateway_key)
    {
    }

    /**
     * Evolve @rule for @epochs
     * @param rule rule number
     * @param epochs number of timesteps to evolve for
     * @return vector of state history over time
     */
    vector<GlobalTransitionOutputType> evolve_rule(size_t rule, size_t epochs)
    {
        cout << "rule:\t" << rule << endl;
        vector<GlobalTransitionOutputType> state_history;

        // Reset state to initial state
        this->_state = this->gateway_key().start_state();
        state_history.push_back(this->_state);

        // Evolve for @epochs and generate state history
        for (size_t epoch = 0; epoch < epochs; epoch++)
            state_history.push_back(global_transition(rule));

        return state_history;
    }

    /**
     * Evolve all rules for @epochs (optionally write images)
     * @param epochs number of timesteps to evolve each rule for
     * @param write_image write PGM to observe visually
     */
    void evolve_all(size_t epochs, bool write_image = false)
    {
        cout << "start state:\t";
        print_vector(this->gateway_key().start_state());
        cout << endl;

        // Enumerate all rules dependent on the total permutations
        for (size_t rule = 0; rule < pow(2, this->gateway_key().total_permutations()); rule++) {
            vector<GlobalTransitionOutputType> state_history = evolve_rule(rule, epochs);

            if (write_image)
                write_pgm(state_history, rule);
        }
    }
};

// TODO: Finish block CA implementation, remove __attribute__((unused))
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class BlockCA1D : public CA1D<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
  private:
    /**
     * Get the block at target cell
     * @param cell target cell
     * @return slice of cells/boundary in partition
     */
    // TODO: assertion on cell mod block size
    vector<CellType> partition(size_t cell)
    {
        // Concatenation of bounded slices from lhs/rhs if at edge, otherwise a standard slice from array
        vector<CellType> partition;

        // TODO: Adjust for biasing
        if (!(this->_gateway_key.partition_size() % 2)) {
            switch (this->_gateway_key.partition_bias()) {
            case PARTITION_BIAS_RHS:
                break;
            case PARTITION_BIAS_CENTER:
                break;
            default:
                break;
            }
        }

        size_t difference = this->_gateway_key.state_size() - cell;
        int bounded_cells = difference - this->_gateway_key.partition_size();

        switch (this->_gateway_key.boundary()) {
        case BOUNDARY_CYCLIC:
            if (bounded_cells < 0) {
                vector<CellType> boundary(this->_state.begin(), this->_state.begin() + abs(bounded_cells));
                vector<CellType> in_bounds(this->_state.begin() + cell, this->_state.end());
                merge_partitions(partition, in_bounds, boundary);
            } else {
                partition =
                    vector<CellType>(this->_state.begin() + cell,
                                     this->_state.begin() + cell + this->_gateway_key.partition_size());
            }
            break;
        case BOUNDARY_ZERO:
            if (bounded_cells < 0) {
                vector<CellType> boundary(abs(bounded_cells), 0);
                vector<CellType> in_bounds(this->_state.begin() + cell, this->_state.end());
                merge_partitions(partition, in_bounds, boundary);
            } else {
                partition =
                    vector<CellType>(this->_state.begin() + cell,
                                     this->_state.begin() + cell + this->_gateway_key.partition_size());
            }
            break;
        default:
            break;
        }

        return partition;
    }

    /**
     * Execute the state transition rule for a block
     * @param cell: the index of the cell in state
     * @param rule: the rule number to be used in the interaction
     * @return: new cell state
     */
    LocalTransitionOutputType local_transition(size_t cell, __attribute__((unused)) size_t rule,
                                               __attribute__((unused)) bool log = true)
    {
        auto partition = this->partition(cell);
        LocalTransitionOutputType new_state;

        cout << "partition:\t";
        print_vector(partition);
        cout << endl;

        // TODO: Trim excess from block added during translation
        switch (this->_gateway_key.interaction()) {
        default:
            cerr << "ERROR: CA1D::local_transition() : unsupported transition" << endl;
            exit(-1);
        }

        return new_state;
    }

    /**
     * Evolve the configured CA for @epochs (global transition function)
     * @param epochs: number of timesteps to run evolution for (entire state update)
     * @param write_image: write PGM files to observe state over time visually
     * @return: ending state
     */
    // TODO: Should this return state history? Or at all?
    GlobalTransitionOutputType global_transition(size_t rule)
    {
        for (size_t block_start = 0; block_start < this->_state.size();
             block_start += this->_gateway_key.partition_size()) {
            LocalTransitionOutputType new_state = local_transition(block_start, rule);

            for (size_t i = block_start; i < PartitionSize; i++)
                this->_state[i] = new_state[i];
        }

        return this->_state;
    }

  public:
    BlockCA1D(GatewayKey<PartitionSize, CellType> gateway_key)
        : CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(gateway_key)
    {
    }

    /**
     * Evolve entire state for one epoch
     * @param rule integer of enumerated rule
     * @return evolved state
     */
    vector<GlobalTransitionOutputType> evolve_rule(size_t rule)
    {
        for (size_t block_start = 0; block_start < this->_state.size();
             block_start += this->_gateway_key.partition_size()) {
            LocalTransitionOutputType new_state = local_transition(block_start, rule);

            // Splice block back into state
            for (size_t i = block_start; i < PartitionSize; i++)
                this->_state[i] = new_state[i];
        }

        return this->_state;
    }

    /**
     * Evolve all enumerated rules
     */
    void evolve_all(size_t epochs, bool write_image)
    {
        for (size_t rule = 0; rule < pow(2, this->_gateway_key.total_permutations()); rule++) {
            cout << "rule:\t" << rule << endl;
            vector<GlobalTransitionOutputType> state_history;

            // Reset state from previous runs
            this->_state = this->_gateway_key.start_state();
            state_history.push_back(this->_state);

            // Evolve
            for (size_t epoch = 0; epoch < epochs; epoch++)
                state_history.push_back(evolve_rule(rule));

            if (write_image)
                write_pgm(state_history, rule);
        }
    }
};

template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class CA2D : public CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
    CA2D(GatewayKey<PartitionSize, CellType> gateway_key)
        : CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(gateway_key)
    {
    }
};

template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class IrreversibleCA2D
    : public CA2D<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
  private:
    /**
     * Get the block at target cell
     * @param cell target cell
     * @return slice of cells/boundary in partition
     */
    // TODO: assertion on cell mod block size
    vector<CellType> partition(__attribute__((unused)) size_t x, __attribute__((unused)) size_t y)
    {
        // Concatenation of bounded slices from lhs/rhs if at edge, otherwise a standard slice from array
        vector<CellType> partition;

        // TODO: Make neighborhood type configurable via gateway key (margolus, moore, etc.)
        switch (this->_gateway_key.boundary()) {
        case BOUNDARY_CYCLIC:
            break;
        case BOUNDARY_ZERO:
            break;
        default:
            break;
        }

        return partition;
    }

    LocalTransitionOutputType local_transition(size_t cell)
    {
        // TODO: Obtain neighbors
        LocalTransitionOutputType partition = this->partition(cell);

        // Logic for Game of Life
        if (this->_state[cell]) {
            // TODO: If not two or three live neighbors, kill
        } else {
            // TODO: If three live neighbors, cell comes alive
        }
    }

    GlobalTransitionOutputType global_transition()
    {
        // TODO: State is 2D, update GatewayKey to accomodate higher dimensions, then break this into x and y
        for (size_t cell = 0; cell < this->state_size(); cell++) {
            this->_state[cell] = local_transition(cell);
        }
    }

  public:
    IrreversibleCA2D(GatewayKey<PartitionSize, CellType> gateway_key)
        : CA2D<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(gateway_key)
    {
    }

    /**
     * Evolve the CA
     */
    void evolve(size_t epochs, bool write_image)
    {
        vector<GlobalTransitionOutputType> state_history;

        // Reset state from previous runs
        this->_state = this->_gateway_key.start_state();
        state_history.push_back(this->_state);

        // Evolve
        for (size_t epoch = 0; epoch < epochs; epoch++)
            state_history.push_back(global_transition());

        // TODO: Update this to create a GIF from bitmaps to observe over time
        // https://github.com/lecram/gifenc
        if (write_image)
            write_pgm(state_history);
    }
};
/// Aliases for well-known/named CAs
namespace Alias1D
{
    typedef IrreversibleCA1D<bool, bool, vector<bool>, 3> ElementaryCA;
    typedef IrreversibleCA2D<bool, bool, vector<vector<bool>>, 8> Life;
}; // namespace Alias1D
