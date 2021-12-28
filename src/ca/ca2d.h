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
            // NOTE: Only collects cells within bounds and excludes
            //       the rest, may need to be expanded for other interactions
            for (int row = y - 1; row <= (int) y + 1; row++) {
                // Exclude out of bound cells
                if (row < 0 || row > (int) this->state().size() - 1) {
                    continue;
                } else {
                    for (int column = x - 1; column <= (int) x + 1; column++) {
                        // Exclude cell itself
                        if ((int) x != column && (int) y != row) {
                            // Exclude out of bound cells
                            if (column < 0 || column > (int) this->_state[0].size() - 1)
                                continue;
                            else
                                partition.push_back(this->_state[row][column]);
                        }
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
                // 2 or 3 neighbors survive, under/overcrowding
                if (living_neighbors != 2 && living_neighbors != 3) {
                    cout << "under/over-crowding" << endl;
                    this->_state[y][x] = false;
                } else {
                    cout << "survive" << endl;
                }
            } else {
                // Reproduction
                if (living_neighbors == 3) {
                    this->_state[y][x] = true;
                    cout << "reproduction" << endl;
                }
            }
            break;
        }
        default:
            cerr << "ERROR: IrreversibleCA2D::local_transition() : unsupported interaction" << endl;
            exit(1);
        }

        return this->_state[y][x];
    }

    // TODO: Optimize (group into blocks, only update those which changed, etc)
    GlobalTransitionOutputType global_transition()
    {
        for (size_t row = 0; row < this->_state.size(); row++) {
            for (size_t column = 0; column < this->_state[0].size(); column++) {
                cout << "cell:\t\t(" << column << "," << row << ")" << endl;
                // Can optimize for early termination by moving neighbors here
                CellType new_state = local_transition(column, row);
                cout << "transition:\t" << this->_state[row][column] << " -> " << new_state << endl;
                this->_state[row][column] = new_state;
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
        GlobalTransitionOutputType last = this->_state;

        if (write_image)
            write_pgm_2d_state(last, 0);

        // Evolve
        for (size_t epoch = 1; epoch < epochs + 1; epoch++) {
            GlobalTransitionOutputType current = global_transition();

            // FIXME: Some still life are oscillators, need more than one timestep
            if (last == current) {
                cout << "converged to still life" << endl;
                break;

                if (write_image)
                    write_pgm_2d_state(current, epoch);
            }

            state_history.push_back(current);
            // TODO: Update this to create a GIF from bitmaps to observe over time
            // https://github.com/lecram/gifenc
            // OR
            // use pyplot
            if (write_image)
                write_pgm_2d_state(current, epoch);
            last = current;
        }
    }
};

namespace Alias2D
{
    typedef IrreversibleCA2D<bool, bool, vector<vector<bool>>, 8> Life;
};
