#include <gtest/gtest.h>
#include <gmock/gmock.h>

// FIXME: Set up include path for tests and reuse compiled objects
#include "../../src/ca/ca1d.h"
#include "../../src/util/states.h"

using namespace ::testing;

// TODO: Cyclic/Zero boundary tests

TEST(CA1D3P, PartitionSizeAndValues)
{
    const size_t PARTITION_SIZE = 3;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>, PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_LHS);
    Alias1D::ElementaryCA ca(gateway_key);

    // TODO: Test boundaries edges instead of slicing entire state
    for (size_t i = 0; i < STATE_SIZE; i++) {
        vector<bool> partition = ca.partition(i);
        // Test the cell is in the correct position of partition
        ASSERT_EQ(partition[1], ca.state()[i]);

        // Test that size is consistent throughout boundary/in-bounds
        ASSERT_EQ(partition.size(), PARTITION_SIZE);
    }
}

// TODO: ----------- Parameterized tests

// TODO: Break into smaller more revealing tests instead of multi-asserting
TEST(CA1D4P, PartitionSizeAndValuesLHS)
{
    const size_t PARTITION_SIZE = 4;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>, PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_LHS);
    IrreversibleCA1D<bool, bool, vector<bool>, PARTITION_SIZE> ca(gateway_key);

    // TODO: Test boundaries edges instead of slicing entire state
    for (size_t i = 0; i < STATE_SIZE; i++) {
        vector<bool> partition = ca.partition(i);
        // Test the cell is in the correct position of partition
        ASSERT_EQ(partition[2], ca.state()[i]);

        // Test that size is consistent throughout boundary/in-bounds
        ASSERT_EQ(partition.size(), PARTITION_SIZE);
    }
}

TEST(CA1D4P, PartitionSizeAndValuesRHS)
{
    const size_t PARTITION_SIZE = 4;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>, PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_RHS);
    IrreversibleCA1D<bool, bool, vector<bool>, PARTITION_SIZE> ca(gateway_key);

    // TODO: Test boundaries edges instead of slicing entire state
    for (size_t i = 0; i < STATE_SIZE; i++) {
        vector<bool> partition = ca.partition(i);
        // Test the cell is in the correct position of partition
        ASSERT_EQ(partition[1], ca.state()[i]);
        // Test that size is consistent throughout boundary/in-bounds
        ASSERT_EQ(partition.size(), PARTITION_SIZE);
    }
}

TEST(CA1D4P, RHSBiasing)
{
    const size_t PARTITION_SIZE = 4;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>, PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_RHS);
    IrreversibleCA1D<bool, bool, vector<bool>, PARTITION_SIZE> ca(gateway_key);
    vector<bool> partition = ca.partition(0);
    ASSERT_THAT(partition, ElementsAre(0, 1, 1, 1));
}

TEST(CA1D4P, LHSBiasing)
{
    const size_t PARTITION_SIZE = 4;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>, PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_LHS);
    IrreversibleCA1D<bool, bool, vector<bool>, PARTITION_SIZE> ca(gateway_key);
    vector<bool> partition = ca.partition(0);
    ASSERT_THAT(partition, ElementsAre(0, 0, 1, 1));
}

// TODO: ----------- END Parameterized tests

TEST(SecondOrder, Reversible)
{
    const size_t STATE_SIZE = 31;
    const size_t EPOCHS = 15;
    vector<bool> start_state = random_1d_start_state(STATE_SIZE);
    vector<bool> prev_state = random_1d_start_state(STATE_SIZE);
    GatewayKey1D<vector<bool>> gateway_key(
        start_state,
        prev_state,
        CA_1D,
        BOUNDARY_CYCLIC,
        PARTITION_BIAS_LHS,
        INTERACTION_SECOND_ORDER_NEIGHBORHOOD_TO_RULE_BIT
    );
    Alias1D::SecondOrderCA ca(gateway_key);
    vector<vector<bool>> forward = ca.evolve_rule(90, EPOCHS);
    vector<vector<bool>> backward = ca.evolve_rule_reverse(90, EPOCHS, forward[forward.size() - 1], forward[forward.size() - 2]);
    vector<vector<bool>> backward_reversed(backward.rbegin(), backward.rend());
    ASSERT_EQ(forward, backward_reversed);
}
