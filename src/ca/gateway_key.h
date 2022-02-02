#pragma once

#include <type_traits>
#include <vector>

using namespace std;

// TODO: Bound interaction types to relevant CA types

// TODO: Specified in template, not exactly used
typedef enum {
    /// Standard 1D cellular automata as displayed in Wolfram's elementary automata
    CA_1D,
    CA_2D
} geometric_shape_t;

typedef enum {
    /// cell = rule_bits[index of partition in config.partition_permutations].
    INTERACTION_NEIGHBORHOOD_TO_RULE_BIT,
    INTERACTION_SECOND_ORDER_NEIGHBORHOOD_TO_RULE_BIT,
    INTERACTION_SECOND_ORDER_NEIGHBORHOOD_TO_RULE_BIT_REVERSE,
    INTERACTION_GAME_OF_LIFE
    // TODO: Create/name rule for block CAs
} interaction_t;

typedef enum {
    /// When boundary is reached, wrap and grab cells from opposite side of partition (also known as periodic)
    BOUNDARY_CYCLIC,
    /// When boundary is reached, use zeroes for needed partition cells
    BOUNDARY_ZERO
} boundary_t;

/// When an even numbered partition is selected, which side of the partition gets the extra cell(s) from
// NOTE: Excludes block CAs, as they do not start from the center in my implementation
typedef enum {
    /// Left
    PARTITION_BIAS_LHS,
    /// Right
    // TODO: Implement for block CAs
    PARTITION_BIAS_RHS,
    // TODO: Used for block CAs where the non-overlapping blocks are wanted to be centered
    PARTITION_BIAS_CENTER
} partition_bias_t;

typedef enum {
    NEIGHBORHOOD_MOORE,
    // TODO
    NEIGHBORHOOD_VON_NEUMANN,
    NEIGHBORHOOD_SMITH,
    NEIGHBORHOOD_COLE,
    NEIGHBORHOOD_MARGOLUS
} neighborhood_t;

// TODO:
typedef enum { UNIFORM, HYBRID } uniformity_t;

// TODO: Create type specific gateway keys. partition permutations are unneeded in higher dimensions, as is
// bias
// TODO: Template start state type
/// Gateway key is synonymous with configuration/encoding of CA.
template<typename StateRepresentation, size_t PartitionSize = 3, typename CellType = bool> class GatewayKey
{
  protected:
    // TODO: Should CA type be removed now that state representation is passed in?
    /// Type of cellular automata (standard, block, second-order, etc.)
    geometric_shape_t _shape;
    /// How to transform cell or neighborhood/block state between timesteps
    interaction_t _interaction;
    /// How to handle needing more cells within the radii at boundary edges
    boundary_t _boundary;
    /// Start state of CA
    StateRepresentation _start_state;
    /// Previous state (used for second order CAs that require two timesteps for initial seed)
    StateRepresentation _prev_state;

    GatewayKey(const StateRepresentation& start_state, geometric_shape_t shape,
               boundary_t boundary = BOUNDARY_ZERO,
               interaction_t interaction = INTERACTION_NEIGHBORHOOD_TO_RULE_BIT)
        : _shape(shape), _interaction(interaction), _boundary(boundary), _start_state(start_state)
    {
    }

    GatewayKey(const StateRepresentation& start_state, const StateRepresentation& prev_state,
               geometric_shape_t shape, boundary_t boundary = BOUNDARY_ZERO,
               interaction_t interaction = INTERACTION_SECOND_ORDER_NEIGHBORHOOD_TO_RULE_BIT)
        : _shape(shape), _interaction(interaction), _boundary(boundary), _start_state(start_state),
          _prev_state(prev_state)
    {
    }

  public:
    geometric_shape_t shape() { return _shape; }

    interaction_t interaction() { return _interaction; }

    boundary_t boundary() { return _boundary; }

    size_t partition_size() { return PartitionSize; }

    StateRepresentation start_state() { return _start_state; }

    StateRepresentation prev_state() { return _prev_state; }

    size_t state_size() { return _start_state.size(); }
};

template<typename StateRepresentation, size_t PartitionSize = 3, typename CellType = bool>
class GatewayKey1D : public GatewayKey<StateRepresentation, PartitionSize, CellType>
{
  private:
    /// Left or right hand preference to pull extra cell from on even-numbered neighborhoods
    partition_bias_t _partition_bias;
    /// All possible permutations of starting neighborhood
    vector<vector<CellType>> _partition_permutations;

    /**
     * @brief Recursive function to generate all permutations of specified partition cells. Used
     * for finding the index of a partition within the exhaustive set in CA1D::*_interaction()
     * and using it as a bit-shift within qualifying rules.
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
    GatewayKey1D(const StateRepresentation& start_state, geometric_shape_t shape,
                 boundary_t boundary = BOUNDARY_ZERO, partition_bias_t partition_bias = PARTITION_BIAS_LHS,
                 interaction_t interaction = INTERACTION_NEIGHBORHOOD_TO_RULE_BIT)
        : GatewayKey<StateRepresentation, PartitionSize, CellType>(start_state, shape, boundary, interaction),
          _partition_bias(partition_bias)
    {
        // Enumerate neighborhood permutations
        vector<CellType> neighborhood(PartitionSize);
        vector<vector<CellType>> permutations;
        generate_partition_states(permutations, neighborhood, 0);
        _partition_permutations = permutations;
    }

    GatewayKey1D(const StateRepresentation& start_state, const StateRepresentation& prev_state,
                 geometric_shape_t shape, boundary_t boundary = BOUNDARY_ZERO,
                 partition_bias_t partition_bias = PARTITION_BIAS_LHS,
                 interaction_t interaction = INTERACTION_NEIGHBORHOOD_TO_RULE_BIT)
        : GatewayKey<StateRepresentation, PartitionSize, CellType>(start_state, prev_state, shape, boundary,
                                                                   interaction),
          _partition_bias(partition_bias)
    {
        // Enumerate neighborhood permutations
        vector<CellType> neighborhood(PartitionSize);
        vector<vector<CellType>> permutations;
        generate_partition_states(permutations, neighborhood, 0);
        _partition_permutations = permutations;
    }

    vector<vector<CellType>> partition_permutations() { return _partition_permutations; }

    size_t total_permutations() { return _partition_permutations.size(); }

    partition_bias_t partition_bias() { return _partition_bias; }
};

// TODO
template<typename StateRepresentation, size_t PartitionSize = 3, typename CellType = bool>
class GatewayKey2D : public GatewayKey<StateRepresentation, PartitionSize, CellType>
{
  private:
    neighborhood_t _neighborhood;

  public:
    GatewayKey2D(const StateRepresentation& start_state, geometric_shape_t shape,
                 boundary_t boundary = BOUNDARY_ZERO, interaction_t interaction = INTERACTION_GAME_OF_LIFE,
                 neighborhood_t neighborhood = NEIGHBORHOOD_MOORE)
        : GatewayKey<StateRepresentation, PartitionSize, CellType>(start_state, shape, boundary, interaction),
          _neighborhood(neighborhood)
    {
    }

    neighborhood_t neighborhood() { return _neighborhood; };
};
