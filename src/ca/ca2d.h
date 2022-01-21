#pragma once

#include <numeric>

#include "../util/special_print.h"
#include "../util/states.h"
#include "ca.h"

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

template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename StateRepresentation = vector<vector<CellType>>, size_t PartitionSize = 8>
class IrreversibleCA2D : public CA2D<CellType, LocalTransitionOutputType, StateRepresentation, PartitionSize>
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
    vector<CellType> partition(int x, int y)
    {
        // TODO: Do any of the other neighborhoods require unflattened types?
        vector<CellType> partition;
        cout << "coordinates:\t" << x << "," << y << endl;
        switch (this->gateway_key().neighborhood()) {
        case NEIGHBORHOOD_MOORE:
            // NOTE: Only collects cells within bounds and excludes
            //       the rest, may need to be expanded for other interactions
            for (int row = y - 1; row <= y + 1; row++) {
                for (int column = x - 1; column <= x + 1; column++) {
                    cout << "N" << row << "," << column << endl;
                    if (row < 0 || row >= (int) this->_state.size()) {
                        partition.push_back(0);
                        cout << "adding zero boundary" << endl;
                        continue;
                    } else {
                        if ((column < 0 || column >= (int) this->_state[0].size())) {
                            partition.push_back(0);
                            cout << "adding zero boundary" << endl;
                        } else if (y == row && x == column) {
                            cout << "skipping cell" << endl;
                            continue;
                        } else {

                            partition.push_back(this->_state[row][column]);
                            cout << "adding " << row << "," << column << endl;
                        }
                    }
                }
            }
            break;
        default:
            cerr << "ERROR: IrreversibleCA2D::partition() : unsupported neighborhood" << endl;
            exit(1);
        }

        cout << "neighbors:\t";
        print_vector(partition);
        cout << endl;
        if (partition.size() != 8) {
            cerr << "Partition size was not 8 @ " << y << "," << x << endl;
            exit(1);
        }
        return partition;
    }

    void local_transition(StateRepresentation& transition, size_t x, size_t y)
    {
        switch (this->gateway_key().interaction()) {
        case INTERACTION_GAME_OF_LIFE: {
            CellType cell = this->_state[y][x];
            vector<CellType> partition = this->partition(x, y);
            int living_neighbors = accumulate(partition.begin(), partition.end(), 0);

            if (/* living */ cell) {
                // 2 or 3 neighbors survive, die from under/overcrowding
                if (living_neighbors <= 1 || living_neighbors >= 4) {
                    cout << "under/over-crowding" << endl;
                    transition[y][x] = false;
                } else {
                    cout << "survival" << endl;
                }
            } else {
                // Reproduction
                if (living_neighbors == 3) {
                    transition[y][x] = true;
                    cout << "reproduction" << endl;
                }
            }
            break;
        }
        default:
            cerr << "ERROR: IrreversibleCA2D::local_transition() : unsupported interaction" << endl;
            exit(1);
        }
    }

    // TODO: Optimize (group into blocks, only update those which changed, etc)
    StateRepresentation global_transition()
    {
        StateRepresentation transition = this->_state;
        for (size_t row = 0; row < this->_state.size(); row++) {
            for (size_t column = 0; column < this->_state[0].size(); column++) {
                cout << "cell:\t\t(" << column << "," << row << ")" << endl;
                // Can optimize for early termination by moving neighbors here
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
     */
    // TODO: Configurable destination
    // TODO: Create directory if it doesn't exist
    // FIXME: Is there a limit on timesteps for oscillators to converge if stuck in oscillation?
    void evolve(size_t epochs, bool write_image)
    {

        // Reset state from previous runs
        this->_state = this->gateway_key().start_state();
        StateRepresentation last = this->_state;

        // Evolve
        for (size_t epoch = 1; epoch < epochs + 1; epoch++) {
            StateRepresentation current = global_transition();

            // FIXME: Some still life are oscillators, need more than one timestep
            if (last == current) {
                cout << "converged to still life" << endl;
                break;
            }

            this->_state = current;

            if (write_image)
                write_pgm(this->_state, epoch);
            last = this->_state;
        }
    }
};

namespace Alias2D
{
    typedef IrreversibleCA2D<bool, bool, vector<vector<bool>>, 8> Life;
};
