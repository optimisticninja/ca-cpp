#include <algorithm>
#include <bitset>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

typedef uint8_t u8;

template<typename T> ostream& operator<<(ostream& out, const vector<T>& v)
{
    out << "{";
    size_t last = v.size() - 1;
    for (size_t i = 0; i < v.size(); ++i) {
        out << std::hex << (int) v[i];
        if (i != last)
            out << ", ";
    }
    out << "}";
    return out;
}

static const bool LINEAR_ADDITIVE_RULE = false;
static const bool NONLINEAR_ADDITIVE_RULE = true;
// linearity = M in paper
static const bool linearity = NONLINEAR_ADDITIVE_RULE;
// Indices for impact coefficients that match the
// variable names in the paper
enum { P, Q, R, S, T };
static vector<vector<bool>> permutation_bases = {
    {0, 0, 0, 0, 1}, {0, 0, 0, 1, 1}, {0, 0, 1, 1, 1}, {0, 1, 1, 1, 1}};

/* Impact coefficients are used to denote connectivity in the
 * additive rule */
vector<bitset<5>> enumerate_impact_coefficients()
{
    vector<bitset<5>> connectivities;
    for (auto& base : permutation_bases) {
        do {
            // Constraint from paper
            if (base[P] & base[T])
                continue;
            bitset<5> x;
            for (auto i = 0; i < base.size(); i++)
                x[i] = base[i];
            connectivities.push_back(x);
        } while (next_permutation(base.begin(), base.end()));
    }

    sort(connectivities.begin(), connectivities.end(),
         [](const auto& lhs, const auto& rhs) { return lhs.to_ulong() < rhs.to_ulong(); });
    return connectivities;
}

int main()
{
    bitset<8> byte(232);
    size_t size = byte.size();
    vector<bitset<5>> connectivities = enumerate_impact_coefficients();
    cout << "Class unspecific connectivities: " << connectivities.size() << endl;
    size_t c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    vector<bitset<5>> class4_rules;
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

    cout << "Class 1 rules: " << c1 << endl;
    cout << "Class 2 rules: " << c2 << endl;
    cout << "Class 3 rules: " << c3 << endl;
    cout << "Class 4 rules: " << c4 << endl;

    vector<u8> random_values;
    for (auto& rule : class4_rules) {
        for (auto j = 0; j < 1000; j++) {
            // Additive rule CA
            random_values.push_back(byte.to_ulong());
            for (auto bit_idx = 0; bit_idx < size; bit_idx++) {
                auto before = byte;
                byte[bit_idx] = linearity ^ rule[P] & byte[(bit_idx - 2) % size] ^
                                rule[Q] & byte[(bit_idx - 1) % size] ^ rule[R] & byte[bit_idx] ^
                                rule[S] & byte[(bit_idx + 1) % size] ^ rule[T] & byte[(bit_idx + 2) % size];
            }
        }
    }

    cout << "Random values: " << random_values << endl;
}
