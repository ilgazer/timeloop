/* Copyright (c) 2019, NVIDIA CORPORATION. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fstream>
#include <iomanip>
#include <utility>

#include "util/accelergy_interface.hpp"

#include "applications/ga-mapper/ga-mapper.hpp"

#include "applications/ga-mapper/mutatable.h"
#include "applications/ga-mapper/random-mapping-generator.hpp"

//--------------------------------------------//
//                Application                 //
//--------------------------------------------//
namespace dse {
  config::CompoundConfigNode find_arch_constraints(config::CompoundConfigNode config)
  {
    // This is primarily a wrapper function written to handle various ways to
    // get to the list of constraints. The only reason there are multiple ways
    // to get to this list is because of backwards compatibility.
    if (config.isList())
    {
      // We're already at the constraints list.
      return config;
    }
    // Constraints can be specified either anonymously as a list, or indirected
    // via a string name.

    std::string name = "";
    if (config.exists("constraints"))
    {
      if (config.lookup("constraints").isList())
      {
        name = "constraints";
      }
      else if (config.lookupValue("constraints", name))
      {
        name = std::string("constraints_") + name;
      }
    }
    else if (config.exists("targets"))
    {
      if (config.lookup("targets").isList())
      {
        name = "targets";
      }
    }

    if (name != "")
    {
      auto constraints = config.lookup(name);
      return constraints;
    }
    return config::CompoundConfigNode(nullptr, YAML::Node(YAML::NodeType::Sequence));
    // else: no constraints specified, nothing to do.
  }

  Application::Application(config::CompoundConfig *config,
                           std::string output_dir_,
                           std::string name) : config_(config), mapspace_(nullptr) {
    auto rootNode = config->getRoot();

    // Version check
    if (rootNode.exists("architecture") && rootNode.lookup("architecture").exists("nodes")) {
      std::cerr << "ERROR: 'nodes' found as a sub-key in the architecture. The 'nodes' key is used by "
          << "the v0.4 timeloopfe front-end format, and will not be recognized by 'timeloop-...' "
          << "commands. Please use the timeloopfe front-end, which is documented at "
          << "https://github.com/Accelergy-Project/timeloop-accelergy-exercises and "
          << "https://timeloop.csail.mit.edu." << std::endl;
      exit(1);
    }

    start = std::chrono::system_clock::now();

    // Model application configuration.
    auto_bypass_on_failure_ = false;
    std::string semi_qualified_prefix = std::move(name);

    if (rootNode.exists("model")) {
      auto model = rootNode.lookup("model");
      model.lookupValue("verbose", verbose_);
      model.lookupValue("auto_bypass_on_failure", auto_bypass_on_failure_);
      model.lookupValue("out_prefix", semi_qualified_prefix);
    }
    std::filesystem::create_directories(output_dir_);
    out_prefix_ = output_dir_ + "/" + semi_qualified_prefix;

    const std::string inds_file_name = out_prefix_ + ".scores.txt";
    inds_file = std::ofstream(inds_file_name);

    // Problem configuration.
    auto problem = rootNode.lookup("problem");
    problem::ParseWorkload(problem, workload_);
    if (verbose_)
      std::cout << "Problem configuration complete." << std::endl;

    // Architecture configuration.
    config::CompoundConfigNode arch;
    if (rootNode.exists("arch")) {
      arch = rootNode.lookup("arch");
    } else if (rootNode.exists("architecture")) {
      arch = rootNode.lookup("architecture");
    }

    bool is_sparse_topology = rootNode.exists("sparse_optimizations");
    arch_specs_ = model::Engine::ParseSpecs(arch, is_sparse_topology);

    if (rootNode.exists("ERT")) {
      auto ert = rootNode.lookup("ERT");
      if (verbose_)
        std::cout << "Found Accelergy ERT (energy reference table), replacing internal energy model." << std::endl;
      arch_specs_.topology.ParseAccelergyERT(ert);
      if (rootNode.exists("ART")) {
        // Nellie: well, if the users have the version of Accelergy that generates ART
        auto art = rootNode.lookup("ART");
        if (verbose_)
          std::cout << "Found Accelergy ART (area reference table), replacing internal area model." << std::endl;
        arch_specs_.topology.ParseAccelergyART(art);
      }
    } else {
#ifdef USE_ACCELERGY
      // Call accelergy ERT with all input files
      if (arch.exists("subtree") || arch.exists("local")) {
        accelergy::invokeAccelergy(config->inFiles, semi_qualified_prefix, output_dir_);
        std::string ertPath = out_prefix_ + ".ERT.yaml";
        auto ertConfig = new config::CompoundConfig(ertPath.c_str());
        auto ert = ertConfig->getRoot().lookup("ERT");
        if (verbose_)
          std::cout << "Generate Accelergy ERT (energy reference table) to replace internal energy model." << std::endl;
        arch_specs_.topology.ParseAccelergyERT(ert);

        std::string artPath = out_prefix_ + ".ART.yaml";
        auto artConfig = new config::CompoundConfig(artPath.c_str());
        auto art = artConfig->getRoot().lookup("ART");
        if (verbose_)
          std::cout << "Generate Accelergy ART (area reference table) to replace internal area model." << std::endl;
        arch_specs_.topology.ParseAccelergyART(art);
      }
#endif
    }

    // Sparse optimizations
    config::CompoundConfigNode sparse_optimizations;
    if (is_sparse_topology)
      sparse_optimizations = rootNode.lookup("sparse_optimizations");
    sparse_optimizations_ = new sparse::SparseOptimizationInfo(
      sparse::ParseAndConstruct(sparse_optimizations, arch_specs_));
    // characterize workload on whether it has metadata
    workload_.SetDefaultDenseTensorFlag(sparse_optimizations_->compression_info.all_ranks_default_dense);

    if (verbose_)
      std::cout << "Sparse optimization configuration complete." << std::endl;

    mapping::load_arch_props(arch_specs_);
    arch_props_ = new ArchProperties(arch_specs_);
    // Architecture constraints.
    config::CompoundConfigNode arch_constraints_root;
    if (arch.exists("constraints"))
      arch_constraints_root = arch.lookup("constraints");
    else if (rootNode.exists("arch_constraints"))
      arch_constraints_root = rootNode.lookup("arch_constraints");
    else if (rootNode.exists("architecture_constraints"))
      arch_constraints_root = rootNode.lookup("architecture_constraints");

    arch_constraints_yaml_ = find_arch_constraints(arch_constraints_root);

    constraints_ = new mapping::Constraints(*arch_props_, workload_);
    constraints_->Parse(arch_constraints_root);

    if (verbose_)
      std::cout << "Architecture configuration complete." << std::endl;

    random_mapping_gen_ = new application::RandomMappingGenerator(config_, arch_constraints_root, &arch_specs_, &workload_);
    // Mapping configuration: expressed as a mapspace or mapping.
    if (!rootNode.lookup("mapping", start_mapping_yaml_))
    {
      start_mapping_yaml_ = random_mapping_gen_->generate_yaml();
    }

    // layout modeling
    std::cout << "Start Parsering Layout" << std::endl;
    config::CompoundConfigNode compound_config_node_layout;
    // ReSharper disable once CppTooWideScope
    bool existing_layout = rootNode.lookup("layout", compound_config_node_layout);

    if (existing_layout) {
      std::map<std::string, std::pair<uint32_t, uint32_t> > externalPortMapping;
      for (std::string &i: arch_specs_.topology.LevelNames())
        externalPortMapping[i] = {
          arch_specs_.topology.GetStorageLevel(i)->num_ports.Get(),
          arch_specs_.topology.GetStorageLevel(i)->num_ports.Get()
        };

      layout_ = layout::ParseAndConstruct(compound_config_node_layout, workload_, externalPortMapping);

      layout_initialized_ = true;

      layout::PrintOverallLayout(layout_);
    } else {
      layout_initialized_ = false;
      std::cout << "No Layout specified, so using bandwidth based modeling" << std::endl;
    }
  }


  Application::~Application() {
    delete mapspace_;
    delete sparse_optimizations_;
  }

  EvalStatus Application::evaluate_ind(IndStats &indStats) {
    if (!indStats.ind.is_valid())
    {
      std::cout << "Failed to satisfy validity" << std::endl;
      return PREFAIL;
    }
    const auto mappingCfg = indStats.ind.to_dict(start_mapping_yaml_);
    indStats.mapping = mapping::ParseAndConstruct(mappingCfg, arch_specs_, workload_);

    bool success = constraints_->SatisfiedBy(&indStats.mapping);
    if (!success) {
      std::cout << "Failed to satisfy constraints. 2\n" << mappingCfg.getYNode() << std::endl;
      return PREFAIL;
    }

    indStats.engine.Spec(arch_specs_);

    if (auto_bypass_on_failure_) {
      const auto pre_eval_status = indStats.engine.
          PreEvaluationCheck(indStats.mapping, workload_, sparse_optimizations_, false);
      for (unsigned level = 1; level < pre_eval_status.size(); level++)
        if (!pre_eval_status[level].success) {
          for (unsigned pvi = 0; pvi < workload_.GetShape()->NumDataSpaces; pvi++)
          {
            // Ugh... mask is offset-by-1 because level 0 is the arithmetic level.
            if (constraints_->BypassStrings()[pvi].at(level - 1) != '1')
            {
              indStats.mapping.datatype_bypass_nest.at(pvi).reset(level - 1);
            }
          }
        }
    }

    auto pre_eval_status = indStats.engine.
        PreEvaluationCheck(indStats.mapping, workload_, sparse_optimizations_, false);
    success = std::accumulate(pre_eval_status.begin(), pre_eval_status.end(), true,
                              [](const bool cur, const model::EvalStatus &status) { return cur && status.success; });
    if (!success) {
      return PREFAIL;
    }

    // Configure the model and evaluate the mapping.
    auto status_per_level = indStats.engine.Evaluate(indStats.mapping, workload_, sparse_optimizations_);
    success = std::accumulate(status_per_level.begin(), status_per_level.end(), true,
                              [](const bool cur, const model::EvalStatus &status) { return cur && status.success; });
    if (success) {
      inds_file << indStats.engine.Energy() << " " << indStats.engine.Cycles() << " ";
      inds_file.flush();
      return SUCCESS;
    }
    return FAIL;
  }

  void Application::save_best_ind_stats() const {
    const std::string config_file_name = out_prefix_ + ".config.yaml";
    const std::string ind_file_name = out_prefix_ + ".best_ind.yaml";
    const std::string stats_file_name = out_prefix_ + ".stats.txt";
    const std::string map_txt_file_name = out_prefix_ + ".map.txt";
    if (!bestIndStats.empty() && bestIndStats[0] != nullptr && bestIndStats[0]->engine.IsEvaluated()) {
      std::ofstream config_file(config_file_name);
      config_file << config_->getYConfig();
      config_file.close();

      std::ofstream ind_file(ind_file_name);
      std::ofstream map_txt_file(map_txt_file_name);
      std::ofstream stats_file(stats_file_name);

      {
        using namespace std::chrono;

        auto stop = system_clock::now();
        auto duration = stop - start;
        auto h = duration_cast<hours>(duration);
        auto m = duration_cast<minutes>(duration % hours(1));
        auto s = duration_cast<seconds>(duration % minutes(1));

        stats_file << "Duration: "
            << std::setfill('0') << std::setw(2) << h.count() << ":"
            << std::setfill('0') << std::setw(2) << m.count() << ":"
            << std::setfill('0') << std::setw(2) << s.count() << std::endl;
      }

      for (auto &indStats: bestIndStats) {
        auto mappingsYaml = YAML::Node();
        auto mappingYaml = YAML::Node();
        mappingYaml["mapping"] = indStats->ind
            .to_dict(config::CompoundConfigNode(nullptr, YAML::Node()))
            .getYNode();
        mappingYaml["stats"]["cycles"] = indStats->engine.Cycles();
        mappingYaml["stats"]["energy"] = indStats->engine.Energy();
        mappingYaml["stats"]["utilisation"] = indStats->engine.Utilization();
        mappingsYaml.push_back(mappingYaml);
        ind_file << mappingsYaml << std::endl;

        indStats->mapping.PrettyPrint(map_txt_file, arch_specs_.topology.StorageLevelNames(),
                                      indStats->engine.GetTopology().UtilizedCapacities(),
                                      indStats->engine.GetTopology().TileSizes());
        map_txt_file << '\n';
        stats_file << indStats->engine << std::endl;
      }
      stats_file.close();
      map_txt_file.close();
      ind_file.close();
      std::cout << std::endl;

      if (!sparse_optimizations_->no_optimization_applied) {
        std::cout << "Summary stats for best mapping found by mapper:" << std::endl;
        std::cout << "  Utilization = " << std::setw(4) << OUT_FLOAT_FORMAT << std::setprecision(2)
            << bestIndStats[0]->engine.Utilization() << " | pJ/Algorithmic-Compute = " << std::setw(8)
            << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << bestIndStats[0]->engine.Energy() /
            static_cast<double>(bestIndStats[0]->engine.GetTopology().AlgorithmicComputes()) << " | pJ/Compute = " <<
            std::setw(8)
            << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << bestIndStats[0]->engine.Energy() /
            static_cast<double>(bestIndStats[0]->engine.GetTopology().ActualComputes()) << std::endl;
      } else {
        std::cout << "Summary stats for best mapping found by mapper:" << std::endl;
        std::cout << "  Energy = " << std::setw(4) << OUT_FLOAT_FORMAT << std::setprecision(2)
            << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << bestIndStats[0]->engine.Energy() << " | pJ/Compute = " <<
            std::setw(8)
            << OUT_FLOAT_FORMAT << PRINTFLOAT_PRECISION << bestIndStats[0]->engine.Energy() /
            static_cast<double>(bestIndStats[0]->engine.GetTopology().ActualComputes()) << std::endl;
      }
    } else {
      std::cout << "MESSAGE: no valid mappings found within search criteria." << std::endl;
    }
  }

  // ---------------
  // Run the mapper.
  // ---------------
  void Application::run() {
    bestIndStats.clear();
    explore_space();

    save_best_ind_stats();
  }
}
