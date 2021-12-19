#include <gtest/gtest.h>
#include <gmock/gmock.h>

// FIXME: Set up include path for tests and reuse compiled objects
#include "../../src/ca/ca1d.h"
#include "../../src/util/states.h"

using namespace ::testing;

// TODO: Cyclic vs Zero boundary tests

TEST(CA1D3P, PartitionSizeAndValues)
{
    const size_t PARTITION_SIZE = 3;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey<PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_LHS);
    CA1D<PARTITION_SIZE, bool> ca(gateway_key);

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
    GatewayKey<PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_LHS);
    CA1D<PARTITION_SIZE, bool> ca(gateway_key);

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
    GatewayKey<PARTITION_SIZE, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_RHS);
    CA1D<PARTITION_SIZE, bool> ca(gateway_key);

    // TODO: Test boundaries edges instead of slicing entire state
    for (size_t i = 0; i < STATE_SIZE; i++) {
        vector<bool> partition = ca.partition(i);
        // Test the cell is in the correct position of partition
        ASSERT_EQ(partition[1], ca.state()[i]);
        // Test that size is consistent throughout boundary/in-bounds
        ASSERT_EQ(partition.size(), PARTITION_SIZE);
    }
}

TEST(CA1DBlock3P, PartitionSizeAndValuesLHS)
{
    const size_t PARTITION_SIZE = 3;
    const size_t STATE_SIZE = 31;
    vector<bool> start_state = wolfram_start_state(STATE_SIZE);
    GatewayKey<PARTITION_SIZE, bool> gateway_key(start_state, CA_1D_BLOCK, BOUNDARY_ZERO, PARTITION_BIAS_RHS);
    CA1D<PARTITION_SIZE, bool> ca(gateway_key);

    // TODO: Test boundaries edges instead of slicing entire state
    for (size_t i = 0; i < STATE_SIZE; i++) {
        vector<bool> partition = ca.partition(i);
        // Test the cell is in the correct position of partition
        ASSERT_EQ(partition[0], ca.state()[i]);
        // Test that size is consistent throughout boundary/in-bounds
        ASSERT_EQ(partition.size(), PARTITION_SIZE);
    }
}

// TODO: CA1DBlock RHS and Center

TEST(CA1D4P, RHSBiasing)
{
    size_t state_size = 31;
    vector<bool> start_state = wolfram_start_state(state_size);
    GatewayKey<4, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_RHS);
    CA1D<4, bool> ca(gateway_key);
    vector<bool> partition = ca.partition(0);
    ASSERT_THAT(partition, ElementsAre(0, 1, 1, 1));
}

TEST(CA1D4P, LHSBiasing)
{
    size_t state_size = 31;
    vector<bool> start_state = wolfram_start_state(state_size);
    GatewayKey<4, bool> gateway_key(start_state, CA_1D, BOUNDARY_ZERO, PARTITION_BIAS_LHS);
    CA1D<4, bool> ca(gateway_key);
    vector<bool> partition = ca.partition(0);
    ASSERT_THAT(partition, ElementsAre(0, 0, 1, 1));
}

// TODO: ----------- END Parameterized tests
