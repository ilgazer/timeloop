//
// Created by ilgaz on 29/10/2025.
//

#include "applications/ga-mapper/nsga_ii.h"

#include <numeric>


namespace dse {
  bool dominates(const IndStats &a, const IndStats &b)
  // do objective values 'a' dominate objective values 'b'?
  // assume positive values and minimisation
  {
    return a.engine.Energy() <= b.engine.Energy() && a.engine.Cycles() <= b.engine.Cycles() &&
           (a.engine.Energy() < b.engine.Energy() || a.engine.Cycles() < b.engine.Cycles());
  }

  std::vector<std::vector<std::unique_ptr<IndStats> > > fast_nd_sort(
    std::vector<std::unique_ptr<IndStats> > &obj_values)
  // fast non-dominated sort
  {
    const int N = static_cast<int>(obj_values.size());

    std::vector<int> pop_ids(N);
    std::iota(pop_ids.begin(), pop_ids.end(), 0);

    std::vector<std::vector<int> > S(N);
    std::vector<int> n(N, 0), rank(N, 0);

    std::vector<std::vector<std::unique_ptr<IndStats> > > fronts(1);
    std::vector<std::vector<int> > front_inds(1);

    for (const int p: pop_ids) {
      S[p].clear(); // S[p] is the set of points dominated by p
      n[p] = 0; // n[p] the number of points that dominate p
      for (int q: pop_ids) {
        if (dominates(*obj_values[p], *obj_values[q])) {
          if (std::find(S[p].begin(), S[p].end(), q) == S[p].end()) //  if (q not in S[p])
            S[p].push_back(q);
        } else if (dominates(*obj_values[q], *obj_values[p])) {
          n[p] += 1;
        }
      }
    }
    for (const int p: pop_ids) {
      // if p is not dominated; add it to the Pareto front
      if (n[p] == 0) {
        // if (p not in front[0])
        if (obj_values[p] != nullptr) {
          obj_values[p]->rank = 0;
          fronts[0].push_back(std::move(obj_values[p]));
          front_inds[0].push_back(p);
        }
      }
    }

    int i = 0;
    while (!fronts[i].empty()) {
      fronts.emplace_back();
      front_inds.emplace_back();
      for (const int p: front_inds[i]) {
        for (int q: S[p]) {
          n[q] -= 1;
          if (n[q] == 0) {
            rank[q] = i + 1;
            if (obj_values[q] != nullptr) {
              fronts[i + 1].push_back(std::move(obj_values[q]));
              front_inds[i + 1].push_back(q);
            }
          }
        }
      }
      i++;
    }

    // delete empty set
    fronts.pop_back();
    return fronts;
  }

  bool compare_energies(const std::unique_ptr<IndStats> &a, const std::unique_ptr<IndStats> &b) {
    return a->engine.Energy() < b->engine.Energy();
  }

  template<typename P, typename Q>
  void _fill_crowding_for_one(std::vector<P> &front, Q getter) {
    std::sort(front.begin(), front.end(), [getter](auto &a, auto &b) { return getter(a) < getter(b); });
    // val = maxval - minval
    const double val = static_cast<double>(getter(front.back())) - static_cast<double>(
                         getter(front[0]));
    const double scale_factor = (val == 0.0) ? 1.0 : val;

    for (std::size_t j = 1; j < front.size() - 1; ++j) {
      const auto next_stat = static_cast<double>(getter(front[j + 1]));
      const auto prev_stat = static_cast<double>(getter(front[j - 1]));
      front[j]->crowding_distance += (next_stat - prev_stat) / scale_factor;
    }
  }

  void fill_crowding_distance(std::vector<std::unique_ptr<IndStats> > &front)
  // calculate crowding distance of the solutions in front
  {
    _fill_crowding_for_one(front, [](const std::unique_ptr<IndStats> &a) { return a->engine.Energy(); });
    _fill_crowding_for_one(front, [](const std::unique_ptr<IndStats> &a) { return a->engine.Cycles(); });
  }
}
