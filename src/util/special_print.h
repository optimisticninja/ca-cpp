#pragma once

#include <iostream>
#include <vector>

using namespace std;

template<typename T> ostream& operator<<(ostream& output, vector<T> const& values)
{
    output << "[ ";
    for (auto const& value : values)
        output << value << " " << endl;
    output << "]";
    return output;
}
