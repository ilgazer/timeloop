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

#pragma once

#include "mutatable.h"
#include "random-mapping-generator.hpp"
#include "mapspaces/mapspace-factory.hpp"
#include "compound-config/compound-config.hpp"
#include "model/sparse-optimization-parser.hpp"
#include "mapping/parser.hpp"
#include "mapping/arch-properties.hpp"
#include "mapping/constraints.hpp"
#include "compound-config/compound-config.hpp"
#include "model/sparse-optimization-parser.hpp"

//--------------------------------------------//
//                Application                 //
//--------------------------------------------//

namespace dse {
  struct IndStats {
    model::Engine engine;
    Mapping mapping;
    Individual ind;
    int rank = 0; //For nsga-ii sorting
    double crowding_distance = 0;
  };

  enum EvalStatus {FAIL, SUCCESS, PREFAIL};

  class Application {
  protected:
    config::CompoundConfig *config_;
    problem::Workload workload_;
    model::Engine::Specs arch_specs_;
    mapspace::MapSpace *mapspace_;
    sparse::SparseOptimizationInfo *sparse_optimizations_;

    std::string out_prefix_ = "timeloop-mapper";
    bool auto_bypass_on_failure_;
    bool verbose_{};
    ArchProperties *arch_props_;
    mapping::Constraints *constraints_{};
    config::CompoundConfigNode start_mapping_yaml_;
    config::CompoundConfigNode arch_constraints_yaml_;
    std::vector<layout::Layout> layout_;
    bool layout_initialized_;

    std::vector<std::unique_ptr<IndStats> > bestIndStats;
    std::ofstream inds_file;
    std::chrono::time_point<std::chrono::system_clock> start;
    application::RandomMappingGenerator * random_mapping_gen_;

  public:
    explicit Application(config::CompoundConfig *config,
                         std::string output_dir_ = ".",
                         std::string name = "timeloop-ga");

    virtual ~Application();

    EvalStatus evaluate_ind(IndStats &indStats);

    virtual void explore_space() = 0;

    void save_best_ind_stats() const;

    // ---------------
    // Run the mapper.
    // ---------------
    void run();

    std::unique_ptr<IndStats> make_valid_random_ind();

  };
}
