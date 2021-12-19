#pragma once

#include <type_traits>
#include <vector>

using namespace std;

typedef enum {
    CA_1D,
    // TODO: below this line
    CA_1D_BLOCK,
    CA_1D_SECOND_ORDER
} ca_t;

typedef enum {
    /// cell = rule_bits[index of partition in config.partition_permutations].
    INTERACTION_NEIGHBORHOOD_TO_RULE_BIT,

    // TODO: below this line
    // FIXME: Dumb rule of my own devising, try to design playground to find interesting rules
    /// cell = rule_bits[index of partition in config.partition_permutations] ^ prev_cell.
    INTERACTION_NEIGHBORHOOD_TO_RULE_BIT_XOR_PREV_CELL
} interaction_t;

typedef enum {
    /// When boundary is reached, wrap and grab cells from opposite side of partition
    BOUNDARY_CYCLIC,
    /// When boundary is reached, use zeroes for needed partition cells
    BOUNDARY_ZERO
} boundary_t;

/// When an even numbered partition is selected, which side of the partition gets the extra cell
// Note: Excludes block CAs, as they do not start from the center in my implementation
typedef enum {
    /// Left
    PARTITION_BIAS_LHS,
    /// Right
    PARTITION_BIAS_RHS
} partition_bias_t;

// TODO: Reorder template parameters for types... values
/// Configurable 1D CAT gateway.
// PartitionSize is synonymous with neighborhood/block dependent on configuration
template<size_t PartitionSize = 3, typename CellType = bool> class GatewayKey
{
  private:
    /// Type of cellular automata (standard, block, second-order, etc.)
    ca_t _ca_type;
    /// How to transform cell or neighborhood/block state between timesteps
    interaction_t _interaction;
    /// How to handle needing more cells within the radii at boundary edges
    boundary_t _boundary;
    /// Left or right hand preference to pull extra cell from on even-numbered neighborhoods
    partition_bias_t _partition_bias;
    /// All possible permutations of starting neighborhood
    vector<vector<CellType>> _partition_permutations;
    /// Start state of CA
    vector<CellType> _start_state;

    /**
     * Recursive function to generate all permutations of specified partition cells
     * @param permutations: vector to store generated permutations in
     * @param partition: starting partition to permute
     * @param cell: cell to operate on
     */
    void generate_partition_states(vector<vector<CellType>>& permutations, vector<CellType>& partition,
                                   size_t cell)
    {
        if (cell == PartitionSize) {
            permutations.push_back(partition);
            return;
        }

        partition[cell] = 0;
        generate_partition_states(permutations, partition, cell + 1);
        partition[cell] = 1;
        generate_partition_states(permutations, partition, cell + 1);
    }

  public:
    GatewayKey(const vector<CellType>& start_state, boundary_t boundary = BOUNDARY_ZERO,
               partition_bias_t partition_bias = PARTITION_BIAS_LHS)
        : // TODO: Make configurable as multiple types are implemented
          _ca_type(CA_1D),
          // TODO: Make configurable as rules are implemented
          _interaction(INTERACTION_NEIGHBORHOOD_TO_RULE_BIT), _boundary(boundary),
          _partition_bias(partition_bias), _start_state(start_state)
    {
        vector<CellType> neighborhood(PartitionSize);
        vector<vector<CellType>> permutations;
        generate_partition_states(permutations, neighborhood, 0);
        _partition_permutations = permutations;
    }

    ca_t ca_type() { return _ca_type; }

    interaction_t interaction() { return _interaction; }

    boundary_t boundary() { return _boundary; }

    partition_bias_t partition_bias() { return _partition_bias; }

    size_t partition_size() { return PartitionSize; }

    vector<vector<CellType>> partition_permutations() { return _partition_permutations; }

    size_t total_permutations() { return _partition_permutations.size(); }

    vector<CellType> start_state() { return _start_state; }

    size_t state_size() { return _start_state.size(); }
};
