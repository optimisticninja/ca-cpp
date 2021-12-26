#pragma once

#include <vector>

#include "gateway_key.h"

using namespace std;

/**
 * @brief abstract cellular automata
 * @tparam LocalTransitionOutputType the output type of the local transition function (cellular primitive or
 * collection for block output)
 * @tparam GlobalTransitionOutputType the output type of the global transition function (type of state
 * representation)
 * @tparam PartitionSize the size of the sliding window used to iterate over state for a transition (including
 * target cell)
 */
template<typename CellType = bool, typename LocalTransitionOutputType = CellType,
         typename GlobalTransitionOutputType = vector<CellType>, size_t PartitionSize = 3>
class CA
{
  protected:
    /// State representation
    GlobalTransitionOutputType _state;

    CA(GlobalTransitionOutputType start_state) : _state(start_state) {}

  public:
    GlobalTransitionOutputType state() { return _state; }
    void state(GlobalTransitionOutputType state) { _state = state; };
};
