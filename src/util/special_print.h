#include <iostream>
#include <vector>

using namespace std;

template<typename T> void print_vector(vector<T> v)
{
    cout << "[ ";
    for (auto c : v) {
        cout << c << " ";
    }
    cout << "]" << endl;
}
