//
// Created by ilgaz on 27/10/2025.
//

#include "applications/ga-mapper/algorithm.h"

#include "applications/ga-mapper/nsga_ii.h"

namespace dse
{
std::string timestamp()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, 1e6 - 1);

  const auto now = std::chrono::system_clock::now();
  const auto in_time_t = std::chrono::system_clock::to_time_t(now);
  std::stringstream ss;
  ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
  ss << '.' << std::setfill('0') << std::setw(6) << dis(gen);
  return ss.str();
}

std::string make_config_path(const config::CompoundConfigNode &dse_config)
{
  std::stringstream ss;
  std::string algorithm;
  assert(dse_config.lookupValue("algorithm", algorithm));
  ss << algorithm;
  std::vector<std::string> mapKeys;
  dse_config.getMapKeys(mapKeys);
  for (auto &key: mapKeys)
  {
    if (key == "algorithm") continue;
    std::string val;
    dse_config.lookupValue(key, val);
    ss << "/" << key << "_" << val;
  }
  return ss.str();
}

std::unique_ptr<Application> parse_config(config::CompoundConfig *config, const std::string &output_dir)
{
  config::CompoundConfigNode dse_config;
  assert(config->getRoot().lookup("dse", dse_config));

  PERM_MUT_PROB = 0.5;
  FACTOR_MUT_PROB = 0.5;
  FACTOR_IMPERFECT_MUT_PROB = 0;
  KEEP_MUT_PROB = 0.3;
  SPLIT_MUT_PROB = 0.5;

  MH_ENABLED = false;
  MH_ROUNDS = INT_MAX;
  MH_SUCCESSFUL_ROUNDS = 2;
  MH_PROPOSAL_STD = 0.4;
  MH_COOLING_SCHEDULE = 1.3;
  MH_CONVERGENCE_ROUNDS = 20000;
  MH_INITIAL_BETA = 4;
  MH_MAX_MUTATION_TRIES = INT_MAX;
  MH_CORRECTION = true;

  dse_config.lookupValue("perm_mut_prob", PERM_MUT_PROB);
  dse_config.lookupValue("factor_mut_prob", FACTOR_MUT_PROB);
  dse_config.lookupValue("factor_imperfect_mut_prob", FACTOR_IMPERFECT_MUT_PROB);
  dse_config.lookupValue("keep_mut_prob", KEEP_MUT_PROB);
  dse_config.lookupValue("split_mut_prob", SPLIT_MUT_PROB);

  dse_config.lookupValue("mh_enabled", MH_ENABLED);
  dse_config.lookupValue("mh_rounds", MH_ROUNDS);
  dse_config.lookupValue("mh_successful_rounds", MH_SUCCESSFUL_ROUNDS);
  dse_config.lookupValue("mh_proposal_std", MH_PROPOSAL_STD);
  dse_config.lookupValue("mh_cooling_schedule", MH_COOLING_SCHEDULE);
  dse_config.lookupValue("mh_convergence_rounds", MH_CONVERGENCE_ROUNDS);
  dse_config.lookupValue("mh_initial_beta", MH_INITIAL_BETA);
  dse_config.lookupValue("mh_max_mutation_tries", MH_MAX_MUTATION_TRIES);
  dse_config.lookupValue("mh_correction", MH_CORRECTION);

  std::string algorithm;
  assert(dse_config.lookupValue("algorithm", algorithm));
  if (algorithm == "hill_climb")
  {
    return std::make_unique<HillClimb>(config, output_dir, "timeloop-ga");
  }
  if (algorithm == "genetic_algorithm")
  {
    return std::make_unique<GeneticAlgorithm>(config, output_dir, "timeloop-ga");
  }
  return nullptr;
}

bool compareInds(const std::unique_ptr<IndStats> &a, const std::unique_ptr<IndStats> &b)
{
  if (!a || !a->engine.IsSpecced())
  {
    return false;
  }
  if (!b || !b->engine.IsSpecced())
  {
    return true;
  }
  return a->engine.Energy() * a->engine.Cycles() < b->engine.Energy() * b->engine.Cycles();
}

std::unique_ptr<IndStats> Application::make_valid_random_ind()
{
  int num_inds = 0;
  // ReSharper disable once CppDFAEndlessLoop
  while (true)
  {
    num_inds += 1;
    auto indStats = std::make_unique<IndStats>(IndStats{
      {}, {}, Individual(workload_, arch_specs_, random_mapping_gen_->generate_yaml(), arch_constraints_yaml_)
    });
    if (!indStats->ind.is_valid())
    {
      continue;
    }
    //Mutation is only useful when we can repair the random mappings.
    if (MH_ENABLED)
    {
      indStats->ind.mutate(workload_);
    }
    indStats->engine.Spec(arch_specs_);
    if (evaluate_ind(*indStats) == SUCCESS)
    {
      std::cout << "tried " << num_inds << std::endl;
      return indStats;
    }
  }
}

HillClimb::HillClimb(config::CompoundConfig *config, std::string output_dir_, std::string name) // NOLINT(*-pro-type-member-init)
  : Application(config, std::move(output_dir_), std::move(name))
{
  config::CompoundConfigNode dse_config;
  assert(config->getRoot().lookup("dse", dse_config));

  dse_config.lookupValue("iters", iters);
  dse_config.lookupValue("auto_bypass_on_failure", auto_bypass_on_failure_);
}

void HillClimb::explore_space()
{
  bestIndStats.push_back(make_valid_random_ind());
  // ReSharper disable once CppDFALoopConditionNotUpdated CppDFAEndlessLoop
  for (int i = 0; i < iters; i++)
  {
    auto indStats = std::make_unique<IndStats>(IndStats{
      {},
      {},
      bestIndStats[0]->ind
    });
    indStats->engine.Spec(arch_specs_);
    if (evaluate_ind(*indStats) != SUCCESS)
    {
      continue;
    }

    if (compareInds(indStats, bestIndStats[0]))
    {
      bestIndStats[0] = std::move(indStats);
    }
  }
}

GeneticAlgorithm::GeneticAlgorithm(config::CompoundConfig *config, std::string output_dir_, // NOLINT(*-pro-type-member-init)
                                   std::string name) : Application(config, std::move(output_dir_), std::move(name))
{
  config::CompoundConfigNode dse_config;
  assert(config->getRoot().lookup("dse", dse_config));

  assert(dse_config.lookupValue("n", n));
  assert(dse_config.lookupValue("gens", gens));
  assert(dse_config.lookupValue("keep", keep));
  first_gen = n;
  dse_config.lookupValue("first_gen", first_gen);
  assert(dse_config.lookupValue("elitism", elitism));
  selection = "truncation";
  dse_config.lookupValue("selection", selection);
  dse_config.lookupValue("auto_bypass_on_failure", auto_bypass_on_failure_);
}

void GeneticAlgorithm::do_selection(
  std::unordered_set<std::unique_ptr<IndStats>, IndStatsPtrHash, IndStatsPtrEq> &thisGenStats)
{
  std::vector<std::unique_ptr<IndStats> > thisGenVector;
  thisGenVector.reserve(thisGenStats.size());
  for (auto it = thisGenStats.begin(); it != thisGenStats.end();)
  {
    // extract() removes the element from the set and returns a node handle
    auto node = thisGenStats.extract(it++);
    thisGenVector.push_back(std::move(node.value()));
  }
  thisGenStats.clear();

  std::random_device rd;
  std::mt19937 gen(rd());

  if (selection == "truncation")
  {
    std::partial_sort(thisGenVector.begin(), thisGenVector.begin() + keep, thisGenVector.end(), compareInds);
    thisGenVector.resize(keep);
    bestIndStats = std::move(thisGenVector);
  }
  else if (selection == "tournament")
  {
    std::uniform_int_distribution thisGenDist(0, static_cast<int>(thisGenVector.size()) - 1);
    bestIndStats.clear();
    for (int j = 0; j < keep; ++j)
    {
      std::unique_ptr<IndStats> *a = nullptr;
      while (a == nullptr || *a == nullptr) a = &thisGenVector[thisGenDist(gen)];
      std::unique_ptr<IndStats> *b = nullptr;
      while (b == nullptr || *b == nullptr) b = &thisGenVector[thisGenDist(gen)];
      if (compareInds(*a, *b))
      {
        bestIndStats.push_back(std::move(*a));
      }
      else
      {
        bestIndStats.push_back(std::move(*b));
      }
    }
  }
  else if (selection == "nsga_ii")
  {
    auto frontiers = fast_nd_sort(thisGenVector);
    std::size_t remaining;
    std::size_t j = 0;
    while ((remaining = static_cast<std::size_t>(keep) - bestIndStats.size()) > 0)
    {
      if (frontiers[j].size() > remaining)
      {
        fill_crowding_distance(frontiers[j]);
        std::partial_sort(
          frontiers[j].begin(),
          frontiers[j].begin() + static_cast<long>(remaining),
          frontiers[j].end(),
          [](auto &a, auto &b) { return a->crowding_distance - b->crowding_distance; }
        );
        std::move(frontiers[j].begin(), frontiers[j].begin() + static_cast<long>(remaining),
                  std::back_inserter(bestIndStats));
      }
      else
      {
        std::move(frontiers[j].begin(), frontiers[j].end(), std::back_inserter(bestIndStats));
      }
      j++;
    }
  }
  else
  {
    std::cout << "Selection option" << selection << "does not exist! Valid options: truncation, tournament, nsga_ii" <<
        "\n";
    exit(1);
  }
}

void GeneticAlgorithm::explore_space()
{
  std::uniform_int_distribution keepDist(0, keep - 1);
  std::random_device rd;
  std::mt19937 gen(rd());

  std::unordered_set<std::unique_ptr<IndStats>, IndStatsPtrHash, IndStatsPtrEq> thisGenStats;

  for (int i = 0; i < first_gen; ++i)
  {
    thisGenStats.insert(make_valid_random_ind());
  }
  do_selection(thisGenStats);

  for (int i = 0; i < gens; i++)
  {
    thisGenStats.clear();
    int all_inds = 0;
    int num_inds = 0;
    inds_file << "\n";

    for (int j = 0; j < n;)
    {
      const auto a = &bestIndStats[keepDist(gen)]->ind;
      const auto b = &bestIndStats[keepDist(gen)]->ind;
      for (auto &ind: a->crossover(*b))
      {
        if (static_cast<int>(thisGenStats.size()) >= n) break;

        auto indStats = std::make_unique<IndStats>(IndStats{{}, {}, std::move(ind)});
        indStats->engine.Spec(arch_specs_);
        indStats->ind.mutate(workload_);
        const auto result = evaluate_ind(*indStats);
        all_inds++;
        if (result == SUCCESS)
        {
          num_inds += 1;
          j++;
          thisGenStats.insert(std::move(indStats));
        }
      }
    }

    if (elitism)
    {
      for (auto &best_ind_stat: bestIndStats)
      {
        thisGenStats.insert(std::move(best_ind_stat));
      }
    }
    bestIndStats.clear();
    do_selection(thisGenStats);
    std::cout << "gen: " << i << ", total: " << num_inds << "/" << all_inds << std::endl;
  }
  inds_file << "\n";
}
}
