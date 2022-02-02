#pragma once

#include <vector>

#include "gateway_key.h"

using namespace std;

/**
 * @brief abstract cellular automata
 * @tparam LocalTransitionOutputType the output type of the local transition function (cellular primitive or
 * collection for block output)
 * @tparam StateRepresentation the output type of the global transition function (type of state
 * representation)
 * @tparam PartitionSize the size of the sliding window used to iterate over state for a transition (including
 * target cell)
 */
template<typename CellType = bool, typename StateRepresentation = vector<CellType>> class CA
{
  protected:
    /// State representation
    StateRepresentation _state;

    CA(StateRepresentation start_state) : _state(start_state) {}

  public:
    StateRepresentation state() { return _state; }
    void state(StateRepresentation state) { _state = state; };
};
