#ifndef SIMULATIONS_CARTESIAN_FF_create_surpass_force_field_HH
#define SIMULATIONS_CARTESIAN_FF_create_surpass_force_field_HH

#include <vector>
#include <iostream>
#include <math.h>
#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <simulations/systems/surpass/SurpassModel.hh>


#include <simulations/forcefields/TotalEnergyByResidue.hh>
#include <simulations/forcefields/surpass/SurpassR12.hh>
#include <simulations/forcefields/surpass/SurpassR13.hh>
#include <simulations/forcefields/surpass/SurpassR14.hh>
#include <simulations/forcefields/surpass/SurpassR15.hh>
#include <simulations/forcefields/surpass/SurpassA13.hh>
#include <simulations/forcefields/surpass/SurpassHydrogenBond.hh>
#include <simulations/forcefields/surpass/SurpassContactEnergy.hh>
#include <simulations/forcefields/surpass/SurpassCentrosymetricEnergy.hh>
#include <simulations/forcefields/surpass/SurpassLocalRepulsionEnergy.hh>
#include <simulations/forcefields/surpass/SurpassHelixStifnessEnergy.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

using core::real;
using simulations::atom_index;
using simulations::residue_index;

static utils::Logger logger("surpass_force_field_factory");

template<class C>
std::shared_ptr<TotalEnergyByResidue> create_surpass_energy(systems::surpass::SurpassModel<C> &system,
    const core::data::sequence::SecondaryStructure_SP ss2, std::istream & input_stream) {

  using namespace simulations::forcefields;

  std::shared_ptr<TotalEnergyByResidue> out = std::make_shared<TotalEnergyByResidue>();

  std::string line;
  std::vector<std::string> tokens;
  std::shared_ptr<ByResidueEnergy> en = nullptr;
  while (std::getline(input_stream, line)) {
    tokens.clear();
    if (line.length() < 5) continue;
    if (line[0] == '#') continue;
    std::replace(line.begin(), line.end(), '\t', ' ');
    utils::split(utils::trim(line), tokens, ' ');
    if (tokens.size() < 2) {
      logger << utils::LogLevel::WARNING << "Skipping a line with less than 2 tokens: " << line << "\n";
      for (std::string &t : tokens)
        logger << "\t\t>" << t << "<\n";
      continue;
    }
    const std::string score_name = tokens.front();
    tokens.erase(tokens.begin());
    const core::real factor = utils::from_string<core::real>(tokens.front());
    tokens.erase(tokens.begin());

    logger << utils::LogLevel::INFO << "Creating " << score_name << " energy function with weight " << factor << "\n";
    logger << utils::LogLevel::FINE << tokens.size() << " arguments passed to its constructor: " <<
    utils::to_string(tokens, " ") << "\n";

    if (score_name.compare("SurpassHydrogenBond") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassHydrogenBond<C>>(system));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassContactEnergy") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassContactEnergy<C>>(system, tokens));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassCentrosymetricEnergy") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassCentrosymetricEnergy < C>>
      (system));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassLocalRepulsionEnergy") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassLocalRepulsionEnergy < C>>
      (system));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassHelixStifnessEnergy") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassHelixStifnessEnergy < C>>
      (system));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassR12") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassR12<C>>(system, ss2, tokens));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassR13") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassR13<C>>(system, ss2, tokens));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassR14") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassR14<C>>(system, ss2, tokens));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassR15") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassR15<C>>(system, ss2, tokens));
      out->add_component(en, factor);
      continue;
    }

    if (score_name.compare("SurpassA13") == 0) {
      en = std::static_pointer_cast<ByResidueEnergy>(std::make_shared<SurpassA13<C>>(system, ss2, tokens));
      out->add_component(en, factor);
      continue;
    }
  }

  return out;

}

template<class C>
std::shared_ptr<TotalEnergyByResidue> create_surpass_energy(systems::surpass::SurpassModel<C> &system,
    const core::data::sequence::SecondaryStructure_SP ss2, const std::string &scfx_cfg) {

  std::stringstream stream(scfx_cfg);
  return create_surpass_energy<C>(system,ss2,stream);
}

} // ~ simulations
} // ~ cartesian
} // ~ ff

#endif
