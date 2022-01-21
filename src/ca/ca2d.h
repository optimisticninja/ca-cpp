#pragma once

#include <numeric>

#include "../util/special_print.h"
#include "../util/states.h"
#include "ca.h"

/**
 * Abstract 2D CA
 * @tparam CellType the type of the cells in the CA
 * @tparam LocalTransitionOutputType the output type of a local update
 * @tparam StateRepresentation the output type of the whole board (global update)
 * @tparam PartitionSize the size of the partition (neighborhood)
 */
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename StateRepresentation = vector<vector<CellType>>, size_t PartitionSize = 8>
class CA2D : public CA<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>
{
  private:
    /// Configuration/encoding of the CA
    GatewayKey2D<StateRepresentation, PartitionSize, CellType> _gateway_key;

  protected:
    CA2D(GatewayKey2D<StateRepresentation, PartitionSize, CellType> gateway_key)
        : CA<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>(
              gateway_key.start_state()),
          _gateway_key(gateway_key)
    {
    }

  public:
    GatewayKey2D<StateRepresentation, PartitionSize, CellType> gateway_key() { return _gateway_key; }
};

/**
 * Irreversible 2D CA implementation
 * @tparam CellType the type of the cells in the CA
 * @tparam LocalTransitionOutputType the output type of a local update
 * @tparam StateRepresentation the output type of the whole board (global update)
 * @tparam PartitionSize the size of the partition (neighborhood)
 */
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename StateRepresentation = vector<vector<CellType>>, size_t PartitionSize = 8>
class IrreversibleCA2D : public CA2D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>
{
  private:
    /**
     * Get the partition (neighborhood) at target cell
     * @param x column of target cell
     * @param y row of target cell
     * @return slice of cells/boundary in partition
     */
    vector<CellType> partition(int x, int y)
    {
        // TODO: Do any of the other neighborhoods require unflattened types?
        vector<CellType> partition;
        switch (this->gateway_key().neighborhood()) {
        case NEIGHBORHOOD_MOORE:
            // NOTE: Only collects cells within bounds and excludes
            //       the rest, may need to be expanded for other interactions
            for (int row = y - 1; row <= y + 1; row++) {
                for (int column = x - 1; column <= x + 1; column++) {
                    // Potential FIXME: Create entire neighborhood instead of skipping
                    // out of bounds for CAs more complex than life
                    if (row < 0 || row >= (int) this->_state.size()) {
                        continue;
                    } else {
                        if ((column < 0 || column >= (int) this->_state[0].size())) {
                            continue;
                        } else if (y == row && x == column) {
                            continue;
                        } else {
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

    /**
     * Local transition of CA (cellular or block level)
     * @param transition reference to the timestep after current that gets updated
     * @param x column of the target cell
     * @param y row of the target cell
     */
    void local_transition(StateRepresentation& transition, size_t x, size_t y)
    {
        switch (this->gateway_key().interaction()) {
        case INTERACTION_GAME_OF_LIFE: {
            CellType cell = this->_state[y][x];
            vector<CellType> partition = this->partition(x, y);
            int living_neighbors = accumulate(partition.begin(), partition.end(), 0);

            if (cell) {
                if (living_neighbors == 2 || living_neighbors == 3) {
                    transition[y][x] = false;
                }
            } else {
                if (living_neighbors == 3) {
                    transition[y][x] = true;
                }
            }
            break;
        }
        default:
            cerr << "ERROR: IrreversibleCA2D::local_transition() : unsupported interaction" << endl;
            exit(1);
        }
    }

    /**
     * Global transition of entire board state
     * @return new state at timestep + 1
     */
    // TODO: Optimize (group into blocks, only update those which changed, etc)
    StateRepresentation global_transition()
    {
        StateRepresentation transition = this->_state;
        for (size_t row = 0; row < this->_state.size(); row++) {
            for (size_t column = 0; column < this->_state[0].size(); column++) {
                cout << "cell:\t\t(" << column << "," << row << ")" << endl;
                local_transition(transition, column, row);
                cout << "transition:\t" << this->_state[row][column] << " -> " << transition[row][column]
                     << endl;
            }
        }

        return transition;
    }

  public:
    IrreversibleCA2D(GatewayKey2D<StateRepresentation, PartitionSize, CellType> gateway_key)
        : CA2D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>(gateway_key)
    {
    }

    /**
     * Evolve the CA
     * @param epochs number of timesteps to limit CA runtime to
     * @param write_image write PGM files of different timesteps
     */
    void evolve(size_t epochs, bool write_image)
    {

        // Reset state from previous runs
        this->_state = this->gateway_key().start_state();
        StateRepresentation last = this->_state;

        // Evolve
        for (size_t epoch = 1; epoch < epochs + 1; epoch++) {
            StateRepresentation current = global_transition();

            // FIXME: Some still life are oscillators, need more than one timestep
            // Create queue that automaticall dequeues when size of 3 has new member
            // enqueued
            if (last == current) {
                cout << "converged to still life" << endl;
                break;
            }

            this->_state = current;

            if (write_image)
                write_pgm_2d_state(this->_state, epoch);

            last = this->_state;
        }
    }
};

namespace Alias2D
{
    typedef IrreversibleCA2D<bool, bool, vector<vector<bool>>, 8> Life;
};
