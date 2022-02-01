#pragma once

#include <vector>

using namespace std;

namespace neighborhood
{

    // NOTE: Only collects cells within bounds and excludes
    //       the rest, may need to be expanded for other interactions
    template<typename CellType>
    vector<CellType> moore_flattened(const vector<vector<CellType>>& state, int x, int y)
    {
        vector<CellType> neighborhood;
        for (int row = y - 1; row <= y + 1; row++) {
            for (int column = x - 1; column <= x + 1; column++) {
                // Potential FIXME: Create entire neighborhood instead of skipping
                // out of bounds for CAs more complex than life
                if (row < 0 || row >= (int) state.size()) {
                    continue;
                } else {
                    if ((column < 0 || column >= (int) state[0].size())) {
                        continue;
                    } else if (y == row && x == column) {
                        continue;
                    } else {
                        neighborhood.push_back(state[row][column]);
                    }
                }
            }
        }
        return neighborhood;
    }
}; // namespace neighborhood
