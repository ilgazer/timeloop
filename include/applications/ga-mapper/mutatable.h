//
// Created by ilgaz on 26/10/2025.
//

#ifndef TIMELOOP_MUTATABLE_H
#define TIMELOOP_MUTATABLE_H

#include <map>
#include <string>
#include <boost/container_hash/hash.hpp>

#include "compound-config/compound-config.hpp"
#include "workload/workload.hpp"
#include <symengine/expression.h>
#include <symengine/lambda_double.h>

#include "model/engine.hpp"

extern double PERM_MUT_PROB;
extern double FACTOR_MUT_PROB;
extern double FACTOR_IMPERFECT_MUT_PROB;
extern double KEEP_MUT_PROB;
extern double SPLIT_MUT_PROB;

extern bool MH_ENABLED;
extern double MH_ROUNDS;
extern double MH_SUCCESSFUL_ROUNDS;
extern double MH_PROPOSAL_STD;
extern double MH_COOLING_SCHEDULE;
extern double MH_CONVERGENCE_ROUNDS;
extern double MH_INITIAL_BETA;
extern double MH_MAX_MUTATION_TRIES;
extern bool MH_CORRECTION;

namespace dse
{
class BypassEntry
{
public:
  std::string target;
  std::string type;
  config::CompoundConfigNode orig_;
  std::map<std::string, bool> keep_by_dataspace;
  std::map<std::string, bool> keep_by_dataspace_constraints;

  BypassEntry(const config::CompoundConfigNode &node, const config::CompoundConfigNode &arch_constraints);

  BypassEntry();

  [[nodiscard]] YAML::Node to_dict() const;

  void mutate();

  friend bool operator==(const BypassEntry &lhs, const BypassEntry &rhs);

  friend bool operator!=(const BypassEntry &lhs, const BypassEntry &rhs);

  friend std::size_t hash_value(const BypassEntry &obj);
};

struct Factor
{
  std::string target;
  std::string type;
  std::string factor_name;
  SymEngine::RCP<const SymEngine::Symbol> symbol;
  std::uint64_t factor;
  std::uint64_t min_factor;
  std::uint64_t max_factor;
  std::uint64_t remainder;

  Factor(const std::string &target, const std::string &type, const std::string &factor_name, const std::uint64_t factor,
         const std::uint64_t min_factor,
         const std::uint64_t max_factor, const std::uint64_t remainder);

  Factor(const std::string &target, const std::string &type, const std::string &factor_name, const std::uint64_t factor,
         const std::uint64_t max_factor,
         const std::uint64_t remainder) : Factor(target, type, factor_name, factor, 1, max_factor, remainder)
  {
  }

  Factor(const std::string &target, const std::string &type, const std::string &factor_name, const std::uint64_t factor,
         const std::uint64_t max_factor) : Factor(target, type, factor_name, factor, max_factor, factor)
  {
  }

  Factor() : Factor("", "", "", 1, 1)
  {
  }

  std::string to_string() const;

  bool is_valid() const;

  friend bool operator==(const Factor &lhs, const Factor &rhs);

  friend bool operator!=(const Factor &lhs, const Factor &rhs);

  friend std::size_t hash_value(const Factor &obj);
};

class LoopEntry
{
public:
  std::string target;
  std::string type;
  std::map<std::string, Factor> factors;
  std::string permutation_prefix;
  std::string permutation;
  std::string permutation_suffix;
  std::optional<std::uint64_t> split;
  bool split_constrained;
  config::CompoundConfigNode orig_;
  std::optional<std::uint64_t> size;
  std::optional<std::uint64_t> fanout;

  void add_permutation_constraints(std::string &buffer);

  LoopEntry(const config::CompoundConfigNode &node, const std::map<std::string, std::uint64_t> &ind_factors,
            const config::CompoundConfigNode &arch_constraints, std::optional<std::uint64_t> size, std::optional<std::uint64_t> fanout);

  LoopEntry(const config::CompoundConfigNode &node, const std::map<std::string, std::uint64_t> &factors,
            const config::CompoundConfigNode &arch_constraints);

  LoopEntry();

  [[nodiscard]] YAML::Node to_dict() const;

  void mutate();

  bool is_valid() const;

  friend bool operator==(const LoopEntry &lhs, const LoopEntry &rhs);

  friend bool operator!=(const LoopEntry &lhs, const LoopEntry &rhs);

  friend std::size_t hash_value(const LoopEntry &obj);
};

class Individual
{
public:
  std::map<unsigned, BypassEntry> bypass_entries;
  std::map<unsigned, LoopEntry> loop_entries;
  std::map<std::string, std::uint64_t> factors;

  Individual(const problem::Workload &workload,
             const model::Engine::Specs &arch_specs,
             const config::CompoundConfigNode &mapping_config,
             const config::CompoundConfigNode &arch_constraint_targets);

  // We have to take the original mapping CompoundConfigNode as input to preserve the reference to the root config
  [[nodiscard]] config::CompoundConfigNode to_dict(config::CompoundConfigNode node) const;

  void fix_factor_remainders();

  void mutate_factors();

  void mutate(problem::Workload &workload);

  SymEngine::Expression make_storage_usage_constraint(problem::Workload &workload,
                                                      std::map<unsigned, LoopEntry>::iterator entry_iter);

  std::tuple<std::vector<Factor *>, std::vector<SymEngine::Expression> > make_constraints(problem::Workload &workload);

  std::pair<std::vector<double>, bool> do_mh_walk(SymEngine::LambdaRealDoubleVisitor &f,
                                                  std::vector<double> cur_assignment);

  void mh_mutate_factors(problem::Workload &workload);

  void crossover_loop_entries(Individual &other);

  void crossover_bypass_entries(Individual &other);

  [[nodiscard]] std::array<Individual, 2> crossover(const Individual &other) const;

  bool is_valid() const;

  friend bool operator==(const Individual &lhs, const Individual &rhs);

  friend bool operator!=(const Individual &lhs, const Individual &rhs);

  friend std::size_t hash_value(const Individual &obj);
};
}

#endif //TIMELOOP_MUTATABLE_H
