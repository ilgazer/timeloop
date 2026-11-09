//
// Created by ilgaz on 27/10/2025.
//

#ifndef TIMELOOP_ALGORITHM_H
#define TIMELOOP_ALGORITHM_H
#include "mutatable.h"
#include "applications/ga-mapper/ga-mapper.hpp"

namespace dse {
  std::unique_ptr<Application> parse_config(config::CompoundConfig *config, const std::string &output_dir);

  class HillClimb : public Application {
    int iters;

  public:
    HillClimb(config::CompoundConfig *config, std::string output_dir_, std::string name);

  private:
    void explore_space() override;
  };

  class GeneticAlgorithm : public Application {
    int n;
    int gens;
    int keep;
    int first_gen;
    bool elitism;
    std::string selection;

  public:
    struct IndStatsPtrEq
    {
      bool operator()(const std::unique_ptr<IndStats> &a, const std::unique_ptr<IndStats> &b) const
      {
        auto eq =  a->ind == b->ind;
        return eq;
      }
    };

    struct IndStatsPtrHash
    {
      std::size_t operator()(const std::unique_ptr<IndStats> &a) const
      {
        if (a == nullptr) return 0;
        auto hv = hash_value(a->ind);
        return hv;
      }
    };

    GeneticAlgorithm(config::CompoundConfig *config, std::string output_dir_, std::string name);

    void do_selection(std::unordered_set<std::unique_ptr<IndStats>, IndStatsPtrHash, IndStatsPtrEq> &thisGenVector);

  private:
    void explore_space() override;
  };
}
#endif //TIMELOOP_ALGORITHM_H
