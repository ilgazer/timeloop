//
// Created by ilgaz on 04/04/2026.
//

#include "applications/ga-mapper/random-mapping-generator.hpp"

#include "mapspaces/mapspace-base.hpp"
#include "mapspaces/mapspace-factory.hpp"

namespace application
{
Mapping RandomMappingGenerator::generate_mapping()
{
  mapspace::ID mapping_id;
  while (true)
  {
    if (!search_.Next(mapping_id))
    {
      std::cout << "We somehow ran out of mappings!" << std::endl;
      exit(1);
    }
    /*
    for (int i = 0; i < static_cast<int>(mapspace::Dimension::Num); ++i)
    {
      std::cout << static_cast<mapspace::Dimension>(i)<< ": "<< mapping_id[i] << std::endl;
    }
    */
    search_.Report(search::Status::Success);
    Mapping mapping(workload_);
    auto construction_status = mapspace_->ConstructMapping(mapping_id, &mapping, false);
    [[maybe_unused]] bool success = std::accumulate(construction_status.begin(), construction_status.end(), true,
                                                    [](bool cur, const mapspace::Status &status)
                                                    {
                                                      return cur && status.success;
                                                    });
    //if (!success) continue;

    return mapping;
  }
}

config::CompoundConfigNode RandomMappingGenerator::generate_yaml()
{
  YAML::Emitter yaml_out;
  yaml_out << YAML::BeginSeq;
  generate_mapping().FormatAsYaml(yaml_out, arch_specs_->topology.StorageLevelNames());
  yaml_out << YAML::EndSeq;
  const auto str = yaml_out.c_str();
  return config::CompoundConfigNode(nullptr, YAML::Load(str), config_root_);
}

RandomMappingGenerator::RandomMappingGenerator(
  config::CompoundConfig *config,
  const config::CompoundConfigNode &arch_constraints,
  model::Engine::Specs *arch_specs,
  problem::Workload *workload
) : mapspace_(mapspace::ParseAndConstruct(config::CompoundConfigNode(nullptr, YAML::Node()), arch_constraints,
                                          *arch_specs, *workload, true)),
    search_(config::CompoundConfigNode(nullptr, YAML::Node()), mapspace_),
    config_root_(config),
    workload_(workload),
    arch_specs_(arch_specs)
{
}
} // application
