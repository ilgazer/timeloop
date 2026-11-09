//
// Created by ilgaz on 04/04/2026.
//

#ifndef TIMELOOP_RANDOMMAPPINGGENERATOR_HPP
#define TIMELOOP_RANDOMMAPPINGGENERATOR_HPP
#include "mapping/mapping.hpp"
#include "mapspaces/mapspace-base.hpp"
#include "search/random.hpp"
#include "search/search.hpp"

namespace application
{
class RandomMappingGenerator
{
  mapspace::MapSpace * mapspace_;
  search::RandomSearch search_;
  config::CompoundConfig * config_root_;
  problem::Workload * workload_;
  model::Engine::Specs * arch_specs_;

public:

  RandomMappingGenerator(config::CompoundConfig *config, const config::CompoundConfigNode &arch_constraints,
                         model::Engine::Specs *arch_specs, problem::Workload *workload);

  Mapping generate_mapping();

  config::CompoundConfigNode generate_yaml();
};
} // application

#endif //TIMELOOP_RANDOMMAPPINGGENERATOR_HPP
