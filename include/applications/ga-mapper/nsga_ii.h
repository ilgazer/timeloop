//
// Created by ilgaz on 29/10/2025.
//

#ifndef TIMELOOP_NSGA_II_H
#define TIMELOOP_NSGA_II_H

#include "model/engine.hpp"
#include "applications/ga-mapper/ga-mapper.hpp"

namespace dse {
  std::vector<std::vector<std::unique_ptr<IndStats> > > fast_nd_sort(
    std::vector<std::unique_ptr<IndStats> > &obj_values);

  void fill_crowding_distance(std::vector<std::unique_ptr<IndStats> > &front);
}
#endif //TIMELOOP_NSGA_II_H
