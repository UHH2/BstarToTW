#pragma once

#include "UHH2/core/include/Event.h"


/** Match object p with genparticle from list of genparticles genparts. A
 * particle is matched to a genparticle if it is closest and is within
 * dR_min of the genparticle Returns NULL if no match was found
 */

template<typename T>
const T* match(const Particle &p, const std::vector<T> &genparts, double dR_min) {
  double closest_dR = DBL_MAX;
  const T* matched = NULL;
  for (auto &genp : genparts)
    {
      double dR = uhh2::deltaR(genp, p);
      if (dR < dR_min && dR < closest_dR)
	{
	  closest_dR = dR;
	  matched = &genp;
	}
    }
  return matched;
}

