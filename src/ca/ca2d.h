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
    // FIXME: Create entire neighborhood instead of skipping out of bounds for rules more complex than
    // life/moore neighborhood, will require making type specific, so look look at "concept" note at top of
    // file. can likely be lifted up to CA type.
    vector<CellType> partition(size_t x, size_t y)
    {
        // TODO: Do any of the other neighborhoods require unflattened types?
        vector<CellType> partition;

        switch (this->gateway_key().neighborhood()) {
        case NEIGHBORHOOD_MOORE:
            for (int i = y - 1; i <= (int) y + 1; i++) {
                if (i < 0 || i > (int) this->state().size() - 1) {
                    continue;
                } else {
                    for (int j = x - 1; j <= (int) x + 1; j++) {
                        if (j < 0 || j > (int) this->_state[0].size() - 1)
                            continue;
                        else
                            partition.push_back(this->_state[i][j]);
                    }
                }
            }

            break;
        default:
            cerr << "ERROR: IrreversibleCA2D::partition() : unsupported neighborhood" << endl;
            exit(1);
        }

        return partition;
    }

    LocalTransitionOutputType local_transition(size_t x, size_t y)
    {
        switch (this->gateway_key().interaction()) {
        case INTERACTION_GAME_OF_LIFE: {
            CellType cell = this->_state[y][x];
            vector<CellType> partition = this->partition(x, y);
            int living_neighbors = 0;

            for (auto cell : partition)
                living_neighbors += cell;

            if (cell) {
                // If not two or three live neighbors, kill
                if (living_neighbors != 2 && living_neighbors != 3)
                    this->_state[y][x] = false;

            } else {
                // If three live neighbors, cell comes alive
                if (living_neighbors == 3)
                    this->_state[y][x] = true;
            }
            break;
        }
        default:
            cerr << "ERROR: CA2D::local_transition() : unsupported interaction" << endl;
            exit(1);
        }

        return this->_state[y][x];
    }

    GlobalTransitionOutputType global_transition()
    {
        for (size_t y = 0; y < this->_state.size(); y++) {
            for (size_t x = 0; x < this->_state[0].size(); x++) {
                cout << "cell:\t\t(" << x << "," << y << ") " << this->_state[y][x];
                CellType new_state = local_transition(x, y);
                cout << " -> " << new_state << endl;
                this->_state[y][x] = new_state;
            }
        }

        return this->_state;
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
        this->_state = this->gateway_key().start_state();
        state_history.push_back(this->_state);

        // Evolve
        for (size_t epoch = 0; epoch < epochs; epoch++) {
            GlobalTransitionOutputType step = global_transition();
            state_history.push_back(step);
            // TODO: Update this to create a GIF from bitmaps to observe over time
            // https://github.com/lecram/gifenc
            if (write_image)
                write_pgm_2d_state(step, epoch);
        }
    }
};

namespace Alias2D
{
    typedef IrreversibleCA2D<bool, bool, vector<vector<bool>>, 8> Life;
};
