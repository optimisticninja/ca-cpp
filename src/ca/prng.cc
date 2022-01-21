#include <algorithm>
#include <bitset>
#include <climits>
#include <cmath>
#include <iostream>
#include <vector>

#include "debug.h"

using namespace std;

/* NOTE:
 * r = radius
 * m = neighborhood = 2r+1
 * l = timesteps
 * n = number of cells
 */

/* NOTE:
 * - A CA characterized by EXOR and/or EXNOR dependence is called an additive CA
 * - If neighborhood dependence is EXOR, it is a non-complemented rule
 * - If neighborhood dependence is EXNOR (inversion of modulo 2 logic), the CA is a complemented rule
 * - If the same rule applies to all cells, the CA is uniform
 * - Conversely, it is called a hybrid CA
 */

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef enum { LINEAR_ADDITIVE_RULE, NONLINEAR_ADDITIVE_RULE } linearity_t;

/// Impact coefficient indices that match the paper
enum { P, Q, R, S, T };

/**
 * @brief 4-neighborhood CA pseudo-random number generator
 * @tparam NumConnectivities max number of bits (neighbors) for generation of impact coefficients
 * (connectivities)
 */
template<size_t NumConnectivities = 5, size_t NumCells = 31> class CAPRNG
{
  private:
    linearity_t _linearity;
    vector<bitset<NumConnectivities>> class4_rules;
    bitset<NumCells> _state;

    /**
     * @brief recursive function to generate cell connectivities (impact coefficients)
     * @param connectivities output vector to store connectivities
     * @param partition starting partition to generate from
     * @param cell index of the target cell in state
     */
    void generate_connectivities(vector<bitset<NumConnectivities>>& connectivities,
                                 bitset<NumConnectivities>& partition, size_t cell)
    {
        if (cell == NumConnectivities) {
            // Constraint from paper, P & T != 1, also skipping "class 0" rules
            if (partition[P] & partition[T] || partition.to_ulong() == 0)
                return;

            // Otherwise add to connectivities
            connectivities.push_back(partition);
            return;
        }

        partition[cell] = 0;
        generate_connectivities(connectivities, partition, cell + 1);
        partition[cell] = 1;
        generate_connectivities(connectivities, partition, cell + 1);
    }

  public:
    CAPRNG(int initial_state = rand(), linearity_t linearity = LINEAR_ADDITIVE_RULE)
        : _linearity(linearity), _state(initial_state)
    {
        vector<bitset<NumConnectivities>> connectivities;

        // Generate impact coefficients (connectivities)
        bitset<NumConnectivities> partition;
        generate_connectivities(connectivities, partition, 0);
        // Sort them
        sort(connectivities.begin(), connectivities.end(),
             [](const auto& lhs, const auto& rhs) { return lhs.to_ulong() < rhs.to_ulong(); });

        size_t c1 = 0, c2 = 0, c3 = 0, c4 = 0;

        for (auto& connectivity : connectivities) {
            switch (connectivity.count()) {
            case 1:
                c1 += 1;
                break;
            case 2:
                c2 += 1;
                break;
            case 3:
                c3 += 1;
                break;
            case 4:
                c4 += 1;
                class4_rules.push_back(connectivity);
                break;
            default:
                break;
            }
        }

        assert(connectivities.size() == 23 && "Number of generated connectivities did not match paper");
        assert(c1 == 5 && "Number of class 1 rules did not match paper");
        assert(c2 == 9 && "Number of class 2 rules did not match paper");
        assert(c3 == 7 && "Number of class 3 rules did not match paper");
        assert(c4 == 2 && "Number of class 4 rules did not match paper");
    }

    bool local_transition(size_t cell)
    {
        // Randomly select a class 4 connectivity
        auto rule = this->class4_rules[rand() % 2];

        // Create periodic indices
        int last_idx = _state.size() - 1;
        int cell_minus_2 = cell - 2;
        cell_minus_2 = cell < 2 ? _state.size() - abs(cell_minus_2) : cell_minus_2;
        int cell_minus_1 = cell - 1;
        cell_minus_1 = cell < 1 ? _state.size() - abs(cell_minus_1) : cell_minus_1;
        int cell_plus_1 = cell + 1;
        cell_plus_1 = cell_plus_1 > last_idx ? cell_plus_1 - last_idx : cell_plus_1;
        int cell_plus_2 = cell + 2;
        cell_plus_2 = cell_plus_2 > last_idx ? cell_plus_2 - last_idx : cell_plus_2;

        return this->_linearity ^ rule[P] & _state[cell_minus_2] ^ rule[Q] & _state[cell_minus_1] ^
               rule[R] & _state[cell] ^ rule[S] & _state[cell_plus_1] ^ rule[T] & _state[cell_plus_2];
    }

    bitset<NumCells> global_transition()
    {
        bitset<NumCells> new_state;

        for (size_t cell = 0; cell < new_state.size(); cell++)
            new_state[cell] = local_transition(cell);

        _state = new_state;
        return new_state;
    }

    bitset<NumCells> state() { return _state; }
};

int main()
{
    srand(time(NULL));

    CAPRNG<5, 31> prng(0xDEADBEEF);
    for (size_t i = 0; i < INT_MAX; i++) {
        cout << prng.global_transition().to_ulong() << endl;
    }
}
