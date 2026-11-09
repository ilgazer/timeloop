//
// Created by ilgaz on 26/10/2025.
//

#include "applications/ga-mapper/mutatable.h"

#include <random>
#include <regex>
#include <sstream>
#include <yaml-cpp/yaml.h>
#include <string>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <symengine/lambda_double.h>
#include <symengine/real_double.h>

#include "mapping/parser.hpp"

double PERM_MUT_PROB;
double FACTOR_MUT_PROB;
double FACTOR_IMPERFECT_MUT_PROB;
double KEEP_MUT_PROB;
double SPLIT_MUT_PROB;

bool MH_ENABLED;
double MH_ROUNDS;
double MH_SUCCESSFUL_ROUNDS;
double MH_PROPOSAL_STD;
double MH_COOLING_SCHEDULE;
double MH_CONVERGENCE_ROUNDS;
double MH_INITIAL_BETA;
double MH_MAX_MUTATION_TRIES;
bool MH_CORRECTION;

namespace dse
{
std::vector<int> get_prime_factors(int num)
{
  assert(num > 0);
  std::vector factors = {1};

  // Handle the factor of 2
  while (num % 2 == 0)
  {
    factors.push_back(2);
    num /= 2;
  }

  // Handle odd factors
  for (int i = 3; i * i <= num; i += 2)
  {
    while (num % i == 0)
    {
      factors.push_back(i);
      num /= i;
    }
  }

  // If the remaining number is a prime greater than 2
  if (num > 1)
  {
    factors.push_back(num);
  }

  return factors;
}

std::map<std::string, Factor> factors_from_string(
  const std::string &target,
  const std::string &type,
  const std::string &factor_string,
  const std::map<std::string, std::uint64_t> &ind_factors,
  [[maybe_unused]] std::optional<std::uint64_t> fanout
)
{
  std::map<std::string, Factor> factors_dict;
  std::stringstream ss(factor_string);
  std::string item;
  std::regex re("([a-zA-Z]+)=?(\\d+)(,(\\d+))?");
  std::smatch match;

  while (ss >> item)
  {
    if (std::regex_match(item, match, re))
    {
      auto k = match[1].str();
      std::uint64_t max_factor = ind_factors.at(k);
      /*
      if (fanout.has_value() && max_factor > fanout.value())
      {
        max_factor = fanout.value();
      }
      */

      if (match[4].matched)
      {
        factors_dict[k] = Factor(
          target,
          type,
          k,
          std::stoi(match[2].str()),
          max_factor,
          std::stoi(match[4].str())
        );
      }
      else
      {
        factors_dict[k] = Factor(
          target,
          type,
          k,
          std::stoi(match[2].str()),
          max_factor
        );
      }
    }
  }
  return factors_dict;
}

void add_constraints_to_factors(
  std::map<std::string, Factor> &factors_dict,
  std::string constr_factor_string,
  const std::map<std::string, std::uint64_t> &ind_factors
)
{
  boost::replace_all(constr_factor_string, ",", " ");
  std::stringstream constr_ss(constr_factor_string);
  std::regex constr_re("([a-zA-Z]+)(>=|<=|==|>|<|=|)(\\d+)");
  std::smatch match;
  std::string item;

  while (constr_ss >> item)
  {
    if (std::regex_match(item, match, constr_re))
    {
      auto k = match[1].str();
      auto op = match[2].str();
      auto val = std::stoi(match[3].str());
      if (val == 0)
      {
        val = ind_factors.at(k);
      }
      if (op == "=" || op == "==" || op == "")
      {
        factors_dict[k].min_factor = val;
        factors_dict[k].max_factor = val;
      }
      else if (op == ">")
      {
        factors_dict[k].min_factor = val + 1;
      }
      else if (op == ">=")
      {
        factors_dict[k].min_factor = val;
      }
      else if (op == "<")
      {
        factors_dict[k].max_factor = val - 1;
      }
      else if (op == "<=")
      {
        factors_dict[k].max_factor = val;
      }
    }
  }
}

std::string factors_to_string(const std::map<std::string, Factor> &factors_dict)
{
  std::stringstream ss;
  bool first = true;
  for (const auto &[k, v]: factors_dict)
  {
    if (!first)
    {
      ss << " ";
    }
    ss << k << "=" << v.to_string();
    first = false;
  }
  return ss.str();
}


BypassEntry::BypassEntry(
  const config::CompoundConfigNode &node,
  const config::CompoundConfigNode &arch_constraints
)
{
  orig_ = config::CompoundConfigNode(node);
  assert(orig_.lookupValue("type", type));
  assert(orig_.lookupValue("target", target));
  std::vector<std::string> keep;
  assert(orig_.lookupArrayValue("keep", keep));
  for (const auto &k: keep)
  {
    keep_by_dataspace[k] = true;
  }
  std::vector<std::string> bypass;
  assert(orig_.lookupArrayValue("bypass", bypass));
  for (const auto &k: bypass)
  {
    keep_by_dataspace[k] = false;
  }

  if (arch_constraints.isList())
  {
    for (int i = 0; i < arch_constraints.getLength(); ++i)
    {
      const auto &arch_node = arch_constraints[i];
      std::string arch_type;
      std::string arch_target;
      assert(arch_node.lookupValue("type", arch_type));
      assert(arch_node.lookupValue("target", arch_target));
      if ((arch_type == "datatype" || arch_type == "bypass" || arch_type == "dataspace") && target == arch_target)
      {
        std::vector<std::string> arch_keep;
        if (node.lookupArrayValue("keep", arch_keep))
        {
          for (const auto &k: arch_keep)
          {
            keep_by_dataspace_constraints[k] = true;
          }
        }

        std::vector<std::string> arch_bypass;
        if (node.lookupArrayValue("bypass", arch_bypass))
        {
          for (const auto &k: arch_bypass)
          {
            keep_by_dataspace_constraints[k] = false;
          }
        }

        std::vector<std::string> mapKeys;
        arch_node.getMapKeys(mapKeys);
        for (const auto &k: mapKeys)
        {
          if (!orig_.exists(k))
          {
            orig_.getYNodeRef()[k] = arch_node.getYNode()[k];
          }
        }
      }
    }
  }
}

BypassEntry::BypassEntry() : orig_(config::CompoundConfigNode(nullptr, YAML::Node()))
{
}

YAML::Node BypassEntry::to_dict() const
{
  YAML::Node root = orig_.getYNode();
  root["target"] = target;
  root["type"] = type;
  root["keep"] = YAML::Node(YAML::NodeType::Sequence);
  root["bypass"] = YAML::Node(YAML::NodeType::Sequence);
  for (const auto &[k, v]: keep_by_dataspace)
  {
    if (v)
    {
      root["keep"].push_back(k);
    }
    else
    {
      root["bypass"].push_back(k);
    }
  }
  return root;
}

void BypassEntry::mutate()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dis(0.0, 1.0);

  for (const auto &[k, v]: keep_by_dataspace)
  {
    //We mutate all keep's that are not constrained with probability KEEP_MUT_PROB
    auto iter = keep_by_dataspace_constraints.find(k);
    if (iter == keep_by_dataspace_constraints.end())
    {
      keep_by_dataspace[k] ^= dis(gen) < KEEP_MUT_PROB;
    }
    else
    {
      keep_by_dataspace[k] = iter->second;
    }
  }
}

Factor::Factor(const std::string &target, const std::string &type, const std::string &factor_name,
               const std::uint64_t factor,
               const std::uint64_t min_factor, const std::uint64_t max_factor,
               const std::uint64_t remainder
) : target(target),
    type(type),
    factor_name(factor_name),
    symbol(SymEngine::symbol(target + "_" + type + "_" + factor_name)),
    factor(factor),
    min_factor(min_factor),
    max_factor(max_factor),
    remainder(remainder)
{
  //assert(max_factor >= factor);
  //assert(factor >= remainder);
  //assert(factor >= min_factor);
}

std::string Factor::to_string() const
{
  return factor == remainder ? std::to_string(factor) : std::to_string(factor) + "," + std::to_string(remainder);
}

bool Factor::is_valid() const
{
  return max_factor >= factor && factor >= min_factor;
}

void LoopEntry::add_permutation_constraints(std::string &buffer)
{
  buffer = buffer.substr(0, buffer.find("#"));

  std::istringstream iss(buffer);
  char token;
  bool separator_seen = false;
  const char separator = '_';
  while (iss >> token)
  {
    if (token == separator)
    {
      separator_seen = true;
      continue;
    }

    if (separator_seen)
    {
      permutation_suffix.push_back(token);
    }
    else
    {
      permutation_prefix.push_back(token);
    }
    permutation = std::regex_replace(permutation, std::regex(std::string(1, token)), "");
  }
}

LoopEntry::LoopEntry(
  const config::CompoundConfigNode &node,
  const std::map<std::string, std::uint64_t> &ind_factors,
  const config::CompoundConfigNode &arch_constraints,
  const std::optional<std::uint64_t> size,
  std::optional<std::uint64_t> fanout
) : size(size), fanout(fanout)
{
  orig_ = config::CompoundConfigNode(node);
  assert(orig_.lookupValue("target", target));
  assert(orig_.lookupValue("type", type));
  std::string factors_str;
  assert(orig_.lookupValue("factors", factors_str));
  factors = factors_from_string(target, type, factors_str, ind_factors, fanout);
  assert(orig_.lookupValue("permutation", permutation));
  unsigned long long _split;
  if (orig_.lookupValue("split", _split))
  {
    split = _split;
  }

  if (arch_constraints.isList())
  {
    for (std::uint64_t i = 0; i < static_cast<std::uint64_t>(arch_constraints.getLength()); ++i)
    {
      const auto &arch_node = arch_constraints[i];
      std::string arch_type;
      std::string arch_target;
      assert(arch_node.lookupValue("type", arch_type));
      assert(arch_node.lookupValue("target", arch_target));
      if (type == arch_type && target == arch_target)
      {
        std::string buffer;
        if (arch_node.lookupValue("permutation", buffer))
        {
          add_permutation_constraints(buffer);
        }

        std::string arch_factors;
        if (arch_node.lookupValue("factors", arch_factors))
        {
          add_constraints_to_factors(factors, arch_factors, ind_factors);
        }

        unsigned long long arch_split;
        if (arch_node.lookupValue("split", arch_split))
        {
          split = arch_split;
          split_constrained = true;
        }

        std::vector<std::string> mapKeys;
        arch_node.getMapKeys(mapKeys);
        for (const auto &k: mapKeys)
        {
          if (!orig_.exists(k))
          {
            orig_.getYNodeRef()[k] = arch_node.getYNode()[k];
          }
        }
      }
    }
  }
}

LoopEntry::LoopEntry(
  const config::CompoundConfigNode &node,
  const std::map<std::string, std::uint64_t> &factors,
  const config::CompoundConfigNode &arch_constraints
) : LoopEntry(node, factors, arch_constraints, std::nullopt, std::nullopt)
{
}

LoopEntry::LoopEntry() : split_constrained(false), orig_(config::CompoundConfigNode(nullptr, YAML::Node()))
{
}

YAML::Node LoopEntry::to_dict() const
{
  YAML::Node root = orig_.getYNode();
  root["target"] = this->target;
  root["type"] = this->type;
  root["factors"] = factors_to_string(this->factors);
  root["permutation"] = this->permutation_prefix + this->permutation + this->permutation_suffix;
  if (this->split.has_value())
  {
    root["split"] = this->split.value();
  }
  return root;
}

void LoopEntry::mutate()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dis(0.0, 1.0);

  while (dis(gen) < PERM_MUT_PROB)
  {
    if (this->permutation.length() < 2) break;
    std::uniform_int_distribution distrib(0, static_cast<int>(this->permutation.length()) - 1);
    int i = 0;
    int j = 0;
    while (i == j)
    {
      i = distrib(gen);
      j = distrib(gen);
    }
    std::swap(this->permutation[i], this->permutation[j]);
  }

  if (this->split.has_value() && !this->split_constrained)
  {
    std::uniform_int_distribution choice_dist(0, 1);
    const std::uint64_t direction = choice_dist(gen) == 0 ? -1 : 1;
    while (dis(gen) < SPLIT_MUT_PROB)
    {
      this->split = this->split.value() + direction;
    }
    this->split = std::clamp<std::uint64_t>(this->split.value(), 0, this->factors.size());
  }
}

bool LoopEntry::is_valid() const
{
  for (const auto &[_, v]: this->factors)
  {
    if (!v.is_valid()) return false;
  }
  return true;
}

Individual::Individual(const problem::Workload &workload,
                       const model::Engine::Specs &arch_specs,
                       const config::CompoundConfigNode &mapping_config,
                       const config::CompoundConfigNode &arch_constraint_targets)
{
  const problem::Shape *shape = workload.GetShape();
  for (auto &[id, name]: shape->FlattenedDimensionIDToName)
  {
    factors[name] = workload.GetFlattenedBound(id);
  }
  assert(mapping_config.isList());

  for (int i = 0; i < mapping_config.getLength(); ++i)
  {
    config::CompoundConfigNode entry = config::CompoundConfigNode(mapping_config[i]);
    std::string type;
    entry.lookupValue("type", type);
    auto tiling_level = mapping::FindTargetTilingLevel(entry, type);

    if (type == "datatype")
    {
      bypass_entries.emplace(
        tiling_level,
        BypassEntry(entry, arch_constraint_targets)
      );
    }
    else if (type == "temporal")
    {
      auto storage = arch_specs.topology.GetStorageLevel(mapping::FindTargetStorageLevel(entry));
      loop_entries.emplace(
        tiling_level,
        LoopEntry(entry, factors, arch_constraint_targets, storage->size.GetOptional(), std::nullopt)
      );
    }
    else if (type == "spatial")
    {
      auto level = mapping::FindTargetStorageLevel(entry);
      if (level >= 1)
      {
        auto instances = arch_specs.topology.GetStorageLevel(level - 1)->instances.GetOptional();
        auto last_instances = arch_specs.topology.GetStorageLevel(level)->instances.GetOptional();
        std::optional<std::uint64_t> fanout = std::nullopt;
        if (instances.has_value() && last_instances.has_value())
        {
          fanout = instances.value() / last_instances.value();
        }
        loop_entries.emplace(
          tiling_level,
          LoopEntry(entry, factors, arch_constraint_targets, std::nullopt, fanout)
        );
      }
    }
    else
    {
      throw std::runtime_error("Unknown type");
    }
  }
}

config::CompoundConfigNode Individual::to_dict(config::CompoundConfigNode node) const
{
  YAML::Node root;
  for (auto &[_, bypass_entry]: bypass_entries)
  {
    root.push_back(bypass_entry.to_dict());
  }
  for (auto &[_, loop_entry]: loop_entries)
  {
    root.push_back(loop_entry.to_dict());
  }
  config::CompoundConfigNode our{nullptr, root};
  node.set(our);
  return node;
}

std::uint64_t ceil_div(const std::uint64_t x, const std::uint64_t y)
{
  return (x + y - 1) / y;
}

void Individual::fix_factor_remainders()
{
  for (const auto &[factor, n_]: factors)
  {
    std::uint64_t n = n_;
    int rest = 1;
    for (auto &[_, loop_entry]: loop_entries)
    {
      rest *= loop_entry.factors[factor].factor;
    }
    bool outermost_loop = true;
    for (auto loop_it = loop_entries.rbegin(); loop_it != loop_entries.rend(); ++loop_it)
    {
      auto &loop = loop_it->second.factors[factor];
      rest /= loop.factor;
      loop.remainder = std::min(ceil_div(n, rest), loop.factor);
      n -= (loop.remainder - 1) * rest;

      if (outermost_loop)
      {
        loop.factor = loop.remainder;
      }
      outermost_loop &= loop.factor == 1;
    }
  }
}

void Individual::mutate_factors()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dis(0.0, 1.0);

  for (const auto &[factor, n]: factors)
  {
    while (dis(gen) < FACTOR_IMPERFECT_MUT_PROB)
    {
      std::uniform_int_distribution entry_dis(0, static_cast<int>(loop_entries.size()) - 1);
      std::bernoulli_distribution bool_dis(0.5);
      const int a_idx = entry_dis(gen);
      int b_idx = entry_dis(gen);
      while (a_idx == b_idx)
      {
        b_idx = entry_dis(gen);
      }
      int p = 1;
      while (dis(gen) < FACTOR_IMPERFECT_MUT_PROB) p++;
      p = bool_dis(gen) ? p : -p;

      LoopEntry &a = loop_entries[a_idx];
      LoopEntry &b = loop_entries[b_idx];

      const auto a_backup = a.factors[factor];
      const auto b_backup = b.factors[factor];

      // We cant go over the constraint bounds for a
      a.factors[factor].factor = std::clamp(a.factors[factor].factor + p, a.factors[factor].min_factor,
                                            a.factors[factor].max_factor);

      b.factors[factor].factor = 1;
      int prod = 1;
      for (auto &[_, loop_entry]: loop_entries)
      {
        prod *= loop_entry.factors[factor].factor;
      }
      b.factors[factor].factor = ceil_div(n, prod);

      // Revert if we exceeded bounds for b
      if (!b.factors[factor].is_valid())
      {
        a.factors[factor] = a_backup;
        b.factors[factor] = b_backup;
      }
    }

    while (dis(gen) < FACTOR_MUT_PROB)
    {
      std::uniform_int_distribution entry_dis(0, static_cast<int>(loop_entries.size()) - 1);
      const int a_idx = entry_dis(gen);
      int b_idx = entry_dis(gen);
      while (a_idx == b_idx)
      {
        b_idx = entry_dis(gen);
      }
      LoopEntry &a = std::next(loop_entries.begin(), a_idx)->second;
      LoopEntry &b = std::next(loop_entries.begin(), b_idx)->second;

      std::vector<int> a_primes = get_prime_factors(a.factors[factor].factor);
      if (a_primes.empty()) continue;

      std::uniform_int_distribution prime_dis(0, static_cast<int>(a_primes.size()) - 1);
      const int chosen = a_primes[prime_dis(gen)];

      const auto a_backup = a.factors[factor];
      const auto b_backup = b.factors[factor];

      a.factors[factor].factor /= chosen;
      b.factors[factor].factor *= chosen;

      // Revert if we exceeded constraint bounds for a or b
      if (!a.factors[factor].is_valid() || !b.factors[factor].is_valid())
      {
        a.factors[factor] = a_backup;
        b.factors[factor] = b_backup;
      }

      std::uint64_t rest = 1;
      for (auto &[_, loop_entry]: loop_entries)
      {
        rest *= loop_entry.factors[factor].factor;
      }
      if (rest < n)
      {
        std::cout << "This should never happen!" << std::endl;
        a.factors[factor] = a_backup;
        b.factors[factor] = b_backup;
      }
    }
  }
  fix_factor_remainders();
}

void Individual::mutate(problem::Workload &workload)
{
  for (auto &[_, entry]: bypass_entries)
  {
    entry.mutate();
  }
  for (auto &[_, entry]: loop_entries)
  {
    entry.mutate();
  }
  if (MH_ENABLED)
  {
    mh_mutate_factors(workload);
  }
  else
  {
    mutate_factors();
  }
}

SymEngine::Expression ge(const SymEngine::Expression &a, const SymEngine::Expression &b)
{
  auto eq = (a + 1) / (b + 1) - 1;
  return SymEngine::pow(SymEngine::min({eq, SymEngine::real_double(0)}), 2);
}

SymEngine::Expression Individual::make_storage_usage_constraint(problem::Workload &workload,
                                                                std::map<unsigned, LoopEntry>::iterator entry_iter)
{
  SymEngine::Expression storage_usage = 0;

  std::unordered_map<std::string, SymEngine::Expression> factors_up_to_here;

  const problem::Shape &shape = *workload.GetShape();
  for (const auto &[_, name]: shape.FactorizedDimensionIDToName)
  {
    factors_up_to_here.try_emplace(name, 1);
  }

  for (auto next_entries_iter = loop_entries.begin(); ; ++next_entries_iter)
  {
    auto next_entry = next_entries_iter->second;
    for (auto &[factor, fobj]: next_entry.factors)
    {
      factors_up_to_here[factor] *= fobj.symbol;
    }
    // To make sure that we don't stop before we hit the current iter.
    if (next_entries_iter == entry_iter)
    {
      break;
    }
  }

  for (auto &[dataspace, keep]: bypass_entries[entry_iter->first].keep_by_dataspace)
  {
    if (keep)
    {
      SymEngine::Expression prod_d2 = 1;
      auto projection = shape.Projections[shape.DataSpaceNameToID.at(dataspace)];
      for (auto &d1: projection)
      {
        SymEngine::Expression sum_d1 = 1;

        for (auto [coeff_id, dim_id]: d1)
        {
          auto factor = shape.FactorizedDimensionIDToName.at(dim_id);
          auto coeff = shape.NumCoefficients == coeff_id ? 1 : workload.GetCoefficient(coeff_id);
          sum_d1 += coeff * (factors_up_to_here.at(factor) - 1);
        }

        prod_d2 *= sum_d1;
      }
      storage_usage += prod_d2;
    }
  }
  return storage_usage;
}

std::tuple<std::vector<Factor *>, std::vector<SymEngine::Expression> > Individual::make_constraints(
  problem::Workload &workload)
{
  std::vector<SymEngine::Expression> constraints;
  std::vector<Factor *> usedFactors;
  SymEngine::map_basic_basic substitutions;
  for (auto &[_, entry]: loop_entries)
  {
    for (auto &[factor_name, factor]: entry.factors)
    {
      if (factor.min_factor == factor.max_factor)
      {
        substitutions[factor.symbol] = SymEngine::real_double(factor.min_factor);
      }
      else
      {
        usedFactors.push_back(&factor);
        constraints.push_back(ge(SymEngine::Expression(factor.symbol), factor.min_factor));
        constraints.push_back(ge(factor.max_factor, SymEngine::Expression(factor.symbol)));
      }
    }
    //auto sub = entry.make_constraints();
    //constraints.insert(constraints.end(), sub.begin(), sub.end());
  }

  for (auto &[key, target]: factors)
  {
    SymEngine::Expression lhs = 1;
    SymEngine::Expression subst_tally = 1;
    std::vector<SymEngine::RCP<const SymEngine::Symbol> > unsubstituted;
    for (auto &[_, entry]: loop_entries)
    {
      lhs *= entry.factors[key].symbol;
      auto result = substitutions.find(entry.factors[key].symbol);
      if (result != substitutions.end())
      {
        subst_tally *= result->second;
      }
      else
      {
        unsubstituted.push_back(entry.factors[key].symbol);
      }
    }
    if (subst_tally.subs({}) == SymEngine::real_double(target))
    {
      for (const auto &symbol: unsubstituted)
      {
        substitutions[symbol] = SymEngine::real_double(1);
      }
    }
    else
    {
      constraints.push_back(ge(lhs, target));
      // To give an upper limit to the search space
      constraints.push_back(ge(1.2 * (target + 1), lhs));
    }
  }

  for (auto entry_iter = loop_entries.begin(); entry_iter != loop_entries.end(); ++entry_iter)
  {
    auto &entry = entry_iter->second;

    if (entry.type == "spatial")
    {
      if (entry.fanout.has_value())
      {
        SymEngine::Expression fanout(1);
        for (const auto &[_, factor]: entry.factors)
        {
          fanout *= factor.symbol;
        }
        constraints.push_back(ge(entry.fanout.value(), fanout));
      }
    }
    else if (entry.type == "temporal")
    {
      if (entry.size.has_value())
      {
        const auto storage_usage = make_storage_usage_constraint(workload, entry_iter);
        constraints.push_back(ge(entry.size.value(), storage_usage));
      }
    }
    // std::cout << entry_iter->second.type << " " << entry_iter->second.target << "\n" << constraints.back() << std::endl;
  }

  // This ensures that we don't bother mutating factors that are predetermined
  for (auto &constraint: constraints)
  {
    constraint = constraint.subs(substitutions);
  }

  return {usedFactors, constraints};
}

// log( Phi(x) ), the log of the standard normal CDF, computed stably
// using erfc rather than erf, since erfc has an accurate asymptotic
// expansion in the tails that erf-based formulas lose to cancellation.
inline double norm_logcdf(double x)
{
  // Phi(x) = 0.5 * erfc(-x / sqrt(2))
  double z = -x / std::sqrt(2.0);
  double e = std::erfc(z);
  if (e > 0.0)
  {
    return std::log(0.5 * e);
  }
  // Deep left tail fallback: asymptotic expansion of erfc for large z
  // erfc(z) ~ exp(-z^2)/(z*sqrt(pi)) * (1 - 1/(2z^2) + ...)
  double z2 = z * z;
  double log_erfc_asym = -z2 - std::log(z * std::sqrt(M_PI))
                         + std::log1p(-1.0 / (2.0 * z2));
  return std::log(0.5) + log_erfc_asym;
}

// log q(y | center) for a single dimension, y already rounded/clamped
double log_q(double y, double center, double sigma)
{
  const double log_center = std::log(center);

  if (y <= 1.0)
  {
    // everything that rounds/clamps to 1 lands in (-inf, 1.5)
    double a = (std::log(1.5) - log_center) / sigma;
    return norm_logcdf(a);
  }

  double a = (std::log(y + 0.5) - log_center) / sigma;
  double b = (std::log(y - 0.5) - log_center) / sigma;

  double log_cdf_a = norm_logcdf(a);
  double log_cdf_b = norm_logcdf(b);

  double diff = log_cdf_b - log_cdf_a;
  return log_cdf_a + std::log1p(-std::exp(diff));
}

// log[ q(cur | prop) / q(prop | cur) ], summed over all dimensions
double hastings_log_correction(const std::vector<double> &cur_assignment,
                               const std::vector<double> &proposal,
                               double sigma)
{
  double total = 0.0;
  for (size_t i = 0; i < cur_assignment.size(); ++i)
  {
    total += log_q(cur_assignment[i], proposal[i], sigma)
        - log_q(proposal[i], cur_assignment[i], sigma);
  }
  return total;
}

std::pair<std::vector<double>, bool> Individual::do_mh_walk(SymEngine::LambdaRealDoubleVisitor &log_prob_fn,
                                                            std::vector<double> cur_assignment)
{
  std::mt19937 rng(std::random_device{}());
  std::normal_distribution<double> normal_dist(0.0, MH_PROPOSAL_STD);
  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

  double beta = MH_INITIAL_BETA;
  int rounds = 0;
  int successful_rounds = 0;
  int convergence_rounds = 0;

  std::vector<double> proposal(cur_assignment.size());
  double log_cur_prob = log_prob_fn.call(cur_assignment);
  while ((rounds < MH_ROUNDS && successful_rounds < MH_SUCCESSFUL_ROUNDS) ||
         (log_cur_prob != 0 && convergence_rounds < MH_CONVERGENCE_ROUNDS))
  {
    for (size_t i = 0; i < cur_assignment.size(); ++i)
    {
      const double raw_assgmt = std::exp(std::log(cur_assignment[i]) + normal_dist(rng));
      proposal[i] = std::max(std::round(raw_assgmt), 1.0);
    }

    const double log_proposal_prob = log_prob_fn.call(proposal);
    const double log_correction = MH_CORRECTION
                                    ? hastings_log_correction(cur_assignment, proposal, MH_PROPOSAL_STD)
                                    : 0.;

    const double log_alpha = beta * (log_proposal_prob - log_cur_prob) + log_correction;

    if (rounds >= MH_ROUNDS || successful_rounds >= MH_SUCCESSFUL_ROUNDS)
    {
      convergence_rounds++;
    }

    if (std::log(uniform_dist(rng)) <= log_alpha)
    {
      cur_assignment = proposal;
      log_cur_prob = log_proposal_prob;

      successful_rounds++;
      beta *= MH_COOLING_SCHEDULE;

      if (rounds >= MH_ROUNDS || successful_rounds >= MH_SUCCESSFUL_ROUNDS)
      {
        beta *= MH_COOLING_SCHEDULE;

        if (log_cur_prob == 0)
        {
          return {cur_assignment, true};
        }
      }
    }
    rounds++;
  }

  return {cur_assignment, false};
}

void Individual::mh_mutate_factors(problem::Workload &workload)
{
  SymEngine::Expression log_prob = 0;
  auto [used_factors, constraints] = make_constraints(workload);
  for (auto constraint: constraints)
  {
    //We multiply it by exp(beta) afterwards since all of our constraints share the same beta.
    log_prob += -constraint;
  }
  SymEngine::LambdaRealDoubleVisitor log_prob_fn;
  SymEngine::vec_basic used_vars;
  std::vector<double> cur_assignment;
  used_vars.reserve(used_factors.size());
  cur_assignment.reserve(used_factors.size());

  SymEngine::map_basic_basic substitutions;

  for (const auto used_factor: used_factors)
  {
    used_vars.push_back(used_factor->symbol);
    cur_assignment.push_back(used_factor->factor);

    substitutions[used_factor->symbol] = SymEngine::real_double(used_factor->factor);
  }
  log_prob_fn.init(used_vars, log_prob);

  /*
  for (auto &constraint: constraints)
  {
    std::cout << constraint << std::endl;
    std::cout << constraint.subs(substitutions) << std::endl;
  }
  */
  int tries = 0;
  bool success = false;
  std::vector<double> new_assignment;
  while (!success)
  {
    if (tries >= MH_MAX_MUTATION_TRIES)
    {
      std::cout << "We failed to mutate after " << tries << " tries!" << std::endl;
      SymEngine::map_basic_basic new_subs;

      for (std::size_t i = 0; i < new_assignment.size(); ++i)
      {
        new_subs[used_vars[i]] = SymEngine::real_double(new_assignment[i]);
        std::cout << SymEngine::Expression(used_vars[i]) << " = " << new_assignment[i] << std::endl;
      }

      for (auto &constraint: constraints)
      {
        std::cout << constraint << std::endl;
        std::cout << constraint.subs(new_subs) << std::endl;
      }
      exit(1);
    }

    std::pair<std::vector<double>, bool> result = do_mh_walk(log_prob_fn, cur_assignment);
    new_assignment = result.first;
    success = result.second;
    tries++;
  }
  std::cout << "Success after try " << tries << std::endl;
  cur_assignment = new_assignment;

  for (std::size_t i = 0; i < cur_assignment.size(); ++i)
  {
    used_factors[i]->factor = static_cast<int>(cur_assignment[i]);
  }

  fix_factor_remainders();
}

void Individual::crossover_loop_entries(Individual &other)
{
  std::map<std::pair<std::string, std::string>, LoopEntry *> a_loop_entries;
  for (auto &[_, e]: this->loop_entries)
  {
    a_loop_entries[{e.target, e.type}] = &e;
  }

  std::map<std::pair<std::string, std::string>, LoopEntry *> b_loop_entries;
  for (auto &[_, e]: other.loop_entries)
  {
    b_loop_entries[{e.target, e.type}] = &e;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dis(0.0, 1.0);

  for (const auto &[factor, _]: factors)
  {
    if (dis(gen) < 0.5)
    {
      for (auto const &[key, val]: a_loop_entries)
      {
        std::swap(a_loop_entries[key]->factors[factor], b_loop_entries[key]->factors[factor]);
      }
    }
  }

  for (auto const &[key, val]: a_loop_entries)
  {
    if (dis(gen) < 0.5)
    {
      std::swap(a_loop_entries[key]->permutation, b_loop_entries[key]->permutation);
    }
  }
}

void Individual::crossover_bypass_entries(Individual &other)
{
  std::map<std::string, BypassEntry *> a_bypass_entries;
  for (auto &[_, e]: this->bypass_entries)
  {
    a_bypass_entries[e.target] = &e;
  }

  std::map<std::string, BypassEntry *> b_bypass_entries;
  for (auto &[_, e]: other.bypass_entries)
  {
    b_bypass_entries[e.target] = &e;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dis(0.0, 1.0);

  for (auto const &[key, val]: a_bypass_entries)
  {
    if (dis(gen) < 0.5)
    {
      std::swap(a_bypass_entries[key]->keep_by_dataspace, b_bypass_entries[key]->keep_by_dataspace);
    }
  }
}

std::array<Individual, 2> Individual::crossover(const Individual &other) const
{
  Individual a = *this;
  Individual b = other;

  a.crossover_loop_entries(b);
  a.crossover_bypass_entries(b);

  return {a, b};
}

bool Individual::is_valid() const
{
  for (const auto &[key, val]: this->loop_entries)
  {
    if (!val.is_valid()) return false;
  }

  for (const auto &[factor, n]: factors)
  {
    std::uint64_t l = 0;

    for (auto loop_it = loop_entries.rbegin(); loop_it != loop_entries.rend(); ++loop_it)
    {
      auto &loop = loop_it->second.factors.at(factor);
      l = l * loop.factor + loop.remainder - 1;
    }
    if (l + 1 != n)
    {
      return false;
    }
  }
  return true;
}

bool operator==(const BypassEntry &lhs, const BypassEntry &rhs)
{
  return lhs.target == rhs.target
         && lhs.type == rhs.type
         && lhs.keep_by_dataspace == rhs.keep_by_dataspace;
}

bool operator!=(const BypassEntry &lhs, const BypassEntry &rhs)
{
  return !(lhs == rhs);
}

std::size_t hash_value(const BypassEntry &obj)
{
  std::size_t seed = 0x7EDC265C;
  boost::hash_combine(seed, obj.target);
  boost::hash_combine(seed, obj.type);
  boost::hash_combine(seed, obj.keep_by_dataspace);
  return seed;
}

bool operator==(const Factor &lhs, const Factor &rhs)
{
  return lhs.factor == rhs.factor;
}

bool operator!=(const Factor &lhs, const Factor &rhs)
{
  return !(lhs == rhs);
}

std::size_t hash_value(const Factor &obj)
{
  std::size_t seed = 0x7BC4A705;
  boost::hash_combine(seed, obj.factor);
  boost::hash_combine(seed, obj.min_factor);
  boost::hash_combine(seed, obj.max_factor);
  boost::hash_combine(seed, obj.remainder);
  return seed;
}

bool operator==(const LoopEntry &lhs, const LoopEntry &rhs)
{
  return lhs.target == rhs.target
         && lhs.type == rhs.type
         && lhs.factors == rhs.factors
         && lhs.permutation == rhs.permutation
         && lhs.split == rhs.split;
}

bool operator!=(const LoopEntry &lhs, const LoopEntry &rhs)
{
  return !(lhs == rhs);
}

std::size_t hash_value(const LoopEntry &obj)
{
  std::size_t seed = 0x31CBB6C4;
  boost::hash_combine(seed, obj.target);
  boost::hash_combine(seed, obj.type);
  boost::hash_combine(seed, obj.factors);
  boost::hash_combine(seed, obj.permutation);
  boost::hash_combine(seed, obj.split);
  return seed;
}

bool operator==(const Individual &lhs, const Individual &rhs)
{
  return lhs.bypass_entries == rhs.bypass_entries
         && lhs.loop_entries == rhs.loop_entries
         && lhs.factors == rhs.factors;
}

bool operator!=(const Individual &lhs, const Individual &rhs)
{
  return !(lhs == rhs);
}

std::size_t hash_value(const Individual &obj)
{
  std::size_t seed = 0x00C6E4B5;
  boost::hash_combine(seed, obj.bypass_entries);
  boost::hash_combine(seed, obj.loop_entries);
  boost::hash_combine(seed, obj.factors);
  return seed;
}
}
