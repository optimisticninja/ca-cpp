#pragma once

#include <vector>

using namespace std;

typedef enum { BIAS_LHS, BIAS_RHS } bias_t;

// template<>
class GatewayKey
{
  private:
    int rule;
    bias_t _bias;
};

template<typename EvolvedType, size_t CASize = 232> class GA
{
  private:
    virtual void permute();
    virtual void rate_fitness();

  public:
    GA();
    void evolve();
};
