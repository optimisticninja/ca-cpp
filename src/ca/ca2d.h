#pragma once

#include "ca.h"

template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<vector<CellType>>, size_t PartitionSize = 8>
class CA2D : public CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
  private:
    /// Configuration/encoding of the CA
    GatewayKey2D<GlobalTransitionOutputType, PartitionSize, CellType> _gateway_key;

  protected:
    CA2D(GatewayKey2D<GlobalTransitionOutputType, PartitionSize, CellType> gateway_key)
        : CA<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>(
              gateway_key.start_state()),
          _gateway_key(gateway_key)

    {
    }

  public:
    GatewayKey2D<GlobalTransitionOutputType, PartitionSize, CellType> gateway_key() { return _gateway_key; }
};

template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<vector<CellType>>, size_t PartitionSize = 8>
class IrreversibleCA2D
    : public CA2D<CellType, LocalTransitionOutputType, GlobalTransitionOutputType, PartitionSize>
{
  private:
    /**
     * Get the block at target cell
     * @param x column of target cell
     * @param y row of target cell
     * @return slice of cells/boundary in partition
     */
    // TODO: assertion on cell mod block size
    vector<CellType> partition(__attribute__((unused)) size_t x, __attribute__((unused)) size_t y)
    {
        // Concatenation of bounded slices from lhs/rhs if at edge, otherwise a standard slice from array
        // TODO: Do any of the other neighborhoods require unflattened types?
        vector<CellType> partition;

        switch (this->gateway_key().neighborhood()) {
        case NEIGHBORHOOD_MOORE:
            // TODO: Obtain nearest neighbors
            break;
        default:
            cerr << "ERROR: IrreversibleCA2D::partition() : unsupported neighborhood" << endl;
            exit(1);
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
    IrreversibleCA2D(GatewayKey2D<GlobalTransitionOutputType, PartitionSize, CellType> gateway_key)
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

namespace Alias2D
{
    typedef IrreversibleCA2D<bool, bool, vector<vector<bool>>, 8> Life;
};
