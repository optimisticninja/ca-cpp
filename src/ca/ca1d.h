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
#include "ca.h"
#include "gateway_key.h"

// TODO: Refactor using the word "partition", apparently there are PCAs, so this term is not a good generic
// TODO: Genericize partition per dimension (use attributes to create blocks or neighborhoods)
// TODO: Use 'concept' from c++20 to restrict to supported types
//     concept Derived = std::is_same<U, Class1>::value || std::is_same<U, Class2>::value;
/// Abstract 1D CA
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename StateRepresentation = vector<CellType>, size_t PartitionSize = 3>
class CA1D : public CA<CellType, StateRepresentation>
{
#ifndef DEBUG
  private:
#endif
    /// Configuration/encoding of the CA
    GatewayKey1D<StateRepresentation, PartitionSize, CellType> _gateway_key;

#ifndef DEBUG
  protected:
#else
  public:
#endif
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

    /**
     * @brief Sliding window of operation across entire state. This is the neighborhood size + 1
     * @param cell index of the target cell in state
     * @return neighborhood including @cell
     */
    vector<CellType> partition(size_t cell)
    {
        int rounded_radius = this->gateway_key().partition_size() / 2;
        // FIXME: LHS is a negative index, leftover from cat-playground port
        int lhs = cell - rounded_radius;
        int rhs = cell + rounded_radius;

        // Concatenation of in/out of bounds slices if at edge, otherwise a standard slice
        vector<CellType> partition;

        // Adjust for biasing
        if (!(this->gateway_key().partition_size() % 2)) {
            switch (this->gateway_key().partition_bias()) {
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
        switch (this->gateway_key().boundary()) {
        case BOUNDARY_CYCLIC:
            if /* Wrap LHS to end of state */ (lhs < 0) {
                vector<CellType> lhs_vec(this->_state.begin() + this->_state.size() + lhs,
                                         this->_state.end());
                vector<CellType> rhs_vec(this->_state.begin(), this->_state.begin() + rhs);
                this->merge_partitions(partition, lhs_vec, rhs_vec);
            } /* Wrap RHS to beginning of state */ else if ((size_t) rhs > this->gateway_key().state_size()) {
                int difference = rhs - this->gateway_key().state_size();
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
            } /* Wrap RHS to beginning of state */ else if ((size_t) rhs > this->gateway_key().state_size()) {
                int difference = rhs - this->gateway_key().state_size();
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

    CA1D(GatewayKey1D<StateRepresentation, PartitionSize, CellType> gateway_key)
        : CA<CellType, StateRepresentation>(gateway_key.start_state()), _gateway_key(gateway_key)
    {
    }

  public:
    GatewayKey1D<StateRepresentation, PartitionSize, CellType> gateway_key() { return _gateway_key; }

    virtual LocalTransitionOutputType local_transition(size_t cell, size_t rule, bool log = false) = 0;
    virtual StateRepresentation global_transition(size_t rule) = 0;
};

// TODO: Create optimized elementary CA for bool only that uses bitsets
// TODO: Condense logging into a single line
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename StateRepresentation = vector<CellType>, size_t PartitionSize = 3>
class IrreversibleCA1D : public CA1D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>
{
#ifdef DEBUG
  public:
#endif

    /**
     * @brief Execute the local state transition rule (partition maps to a permutation index used for getting
     * that bit in the rule)
     * @param cell the index of the cell in state
     * @param rule the rule number to be used in the interaction
     * @param log log single-cell state changes
     * @return: new cell state
     */
    LocalTransitionOutputType local_transition(size_t cell, size_t rule, bool log = false)
    {
        auto partition = this->partition(cell);
        cout << "cell:\t\t" << cell << endl;
        cout << "partition:\t" << partition << endl;

        switch (this->gateway_key().interaction()) {
        // Rule Wolfram uses
        case INTERACTION_NEIGHBORHOOD_TO_RULE_BIT: {
            // get bit offset of partition in this->_gateway_key.partition_permutations
            auto permutations = this->gateway_key().partition_permutations();
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
     * @brief Evolve the full state from one timestep to another
     * @param rule rule number to use in local transition
     * @return ending state
     */
    StateRepresentation global_transition(size_t rule)
    {
        cout << "start state:\t" << this->gateway_key().start_state() << endl;
        StateRepresentation new_state(this->_state.begin(), this->_state.end());

        for (size_t cell = 0; cell < new_state.size(); cell++)
            new_state[cell] = local_transition(cell, rule);

        // State transition
        this->state(new_state);
        cout << "state: " << this->_state << endl;
        return this->_state;
    }

#ifndef DEBUG
  public:
#endif
    IrreversibleCA1D(GatewayKey1D<StateRepresentation, PartitionSize, CellType> gateway_key)
        : CA1D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>(gateway_key)
    {
    }

    /**
     * @brief Evolve @rule for @epochs
     * @param rule rule number
     * @param epochs number of timesteps to evolve for
     * @return vector of state history over time
     */
    vector<StateRepresentation> evolve_rule(size_t rule, size_t epochs)
    {
        cout << "rule:\t" << rule << endl;
        vector<StateRepresentation> state_history;

        // Reset state to initial state
        this->_state = this->gateway_key().start_state();
        state_history.push_back(this->_state);

        // Evolve for @epochs and generate state history
        for (size_t epoch = 0; epoch < epochs; epoch++)
            state_history.push_back(global_transition(rule));

        return state_history;
    }

    /**
     * @brief Evolve all rules for @epochs (optionally write images)
     * @param epochs number of timesteps to evolve each rule for
     * @param write_image write PGM to observe visually
     */
    void evolve_all(size_t epochs, bool write_image = false)
    {
        cout << "start state:\t" << this->gateway_key().start_state() << endl;

        // Enumerate all rules dependent on the total permutations
        for (size_t rule = 0; rule < pow(2, this->gateway_key().total_permutations()); rule++) {
            vector<StateRepresentation> state_history = evolve_rule(rule, epochs);

            if (write_image)
                write_pgm(state_history, rule);
        }
    }
};

// TODO: Create optimized elementary CA for bool only that uses bitsets
// TODO: Condense logging into a single line
/// Reversible 1D CA implementation (second-order and eventually block)
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename StateRepresentation = vector<CellType>, size_t PartitionSize = 3>
class ReversibleCA1D : public CA1D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>
{
#ifdef DEBUG
  public:
#else
  private:
#endif

    vector<CellType> prev_state;

    /**
     * @brief Execute the local state transition rule (partition maps to a permutation index used for getting
     * that bit in the rule)
     * @param cell the index of the cell in state
     * @param rule the rule number to be used in the interaction
     * @param log log single-cell state changes
     * @return: new cell state
     */
    LocalTransitionOutputType local_transition(size_t cell, size_t rule, bool log = false)
    {
        auto partition = this->partition(cell);
        cout << "cell:\t\t" << cell << endl;
        cout << "partition:\t" << partition << endl;

        switch (this->gateway_key().interaction()) {
        case INTERACTION_SECOND_ORDER_NEIGHBORHOOD_TO_RULE_BIT: {
            // get bit offset of partition in this->_gateway_key.partition_permutations
            auto permutations = this->gateway_key().partition_permutations();
            auto it = find(permutations.begin(), permutations.end(), partition);
            auto bit_offset = distance(permutations.begin(), it);
            auto first_order_state = (rule >> bit_offset) & 1;
            LocalTransitionOutputType new_state = prev_state[cell] != first_order_state ? true : false;

            if (log)
                cout << "update:\t\t" << this->_state[cell] << " -> " << first_order_state << endl;
            return new_state;
        }
        default:
            cerr << "ERROR: ReversibleCA1D::local_transition(): unimplemented transition" << endl;
            exit(1);
            break;
        }
    }

    /**
     * @brief Evolve the full state from one timestep to another
     * @param rule rule number to use in local transition
     * @return ending state
     */
    StateRepresentation global_transition(size_t rule)
    {
        cout << "start state 1:\t" << prev_state << endl;
        cout << "start state 2:\t" << this->_state << endl;

        StateRepresentation new_state(this->_state.begin(), this->_state.end());

        for (size_t cell = 0; cell < new_state.size(); cell++)
            new_state[cell] = local_transition(cell, rule);

        // State transition
        prev_state = this->_state;
        this->state(new_state);
        cout << "state: " << this->_state << endl;
        return this->_state;
    }

#ifndef DEBUG
  public:
#endif
    ReversibleCA1D(GatewayKey1D<StateRepresentation, PartitionSize, CellType> gateway_key)
        : CA1D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>(gateway_key),
          prev_state(gateway_key.prev_state())
    {
    }

    /**
     * @brief Evolve @rule for @epochs
     * @param rule rule number
     * @param epochs number of timesteps to evolve for
     * @return vector of state history over time
     */
    vector<StateRepresentation> evolve_rule(size_t rule, size_t epochs)
    {
        cout << "rule:\t" << rule << endl;
        vector<StateRepresentation> state_history;
        state_history.push_back(prev_state);
        state_history.push_back(this->_state);

        // Evolve for @epochs and generate state history
        for (size_t epoch = 0; epoch < epochs; epoch++)
            state_history.push_back(global_transition(rule));

        return state_history;
    }

    // TODO: Feels ugly passing as params, should we just add setters for state/prev?
    /**
     * @brief Evolve @rule for @epochs
     * @param rule rule number
     * @param epochs number of timesteps to evolve for
     * @param last last state of previously evolved CA
     * @param before_last second to last state of previously evolved CA
     * @return vector of state history over time
     */
    vector<StateRepresentation> evolve_rule_reverse(size_t rule, size_t epochs, const vector<bool>& last,
                                                    const vector<bool>& before_last)
    {
        this->_state = before_last;
        this->prev_state = last;
        return evolve_rule(rule, epochs);
    }

    /**
     * @brief Evolve all rules for @epochs (optionally write images)
     * @param epochs number of timesteps to evolve each rule for
     * @param write_image write PGM to observe visually
     */
    void evolve_all(size_t epochs, bool write_image = false)
    {
        cout << "start state:\t" << this->gateway_key().start_state() << endl;

        // Enumerate all rules dependent on the total permutations
        for (size_t rule = 0; rule < pow(2, this->gateway_key().total_permutations()); rule++) {
            vector<StateRepresentation> state_history = evolve_rule(rule, epochs);

            if (write_image)
                write_pgm(state_history, rule, "second-order");
        }
    }
};

/// Aliases for well-known/named CAs
namespace Alias1D
{
    typedef IrreversibleCA1D<bool, bool, vector<bool>, 3> ElementaryCA;
    typedef ReversibleCA1D<bool, bool, vector<bool>, 3> SecondOrderCA;
}; // namespace Alias1D
