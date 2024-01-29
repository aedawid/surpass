#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <map>

#include <core/SURPASSenvironment.hh>
#include <core/calc/numeric/interpolators.hh>
#include <core/calc/numeric/Interpolate1D.hh>
#include <utils/io_utils.hh>
#include <utils/options/input_options.hh>

#include <simulations/forcefields/mf/BoundedMFComponent.hh>
#include <simulations/forcefields/mf/MeanFieldDistributions.hh>

namespace simulations {
namespace forcefields {
namespace mf {

using namespace core::calc::numeric;

std::shared_ptr<MeanFieldDistributions> load_1D_distributions(const std::string & ff_file, const core::real pseudocounts_fraction) {

  utils::Logger logger("load_1D_distributions");
  logger << utils::LogLevel::FILE << "Reading ff file: " << ff_file << "\n";
  std::ifstream in(core::SURPASSenvironment::from_file_or_db(ff_file));
  std::string line, key;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<std::string> header_tokens;
  // ------- The first line is always a file header  - just a comment that holds energy term name
  std::getline(in, line);
  std::string name(line.substr(2));
  core::index2 cnt = 0; // double-check : count the distributions and compare with #keys

  CatmullRomInterpolator<double> cri; // --- interpolator algorithm used by spline interpolated functions
  std::shared_ptr<MeanFieldDistributions> mf_sp = std::make_shared<MeanFieldDistributions>();
  mf_sp->name(name);
  // ------- The second line provides the arguments for an interpolated function
  std::getline(in, line);
  utils::split(line, x);
  core::real left_bound = x.front();
  core::real right_bound = x.back();
  while (std::getline(in, line)) {
    if (line.length() < 2) continue;
    if (line[0] == '#') {
      header_tokens.clear();
      utils::split(line.substr(2),header_tokens);
      key = header_tokens[0];
      if(header_tokens.size()>2) {
        left_bound = utils::from_string<core::real>(header_tokens[1]);
        right_bound = utils::from_string<core::real>(header_tokens[2]);
        std::cerr << key<<" "<<left_bound<<" "<<right_bound<<"\n";
      }
    }
    else {
      y.clear();
      utils::split(line, y);
      if (pseudocounts_fraction > 0) {
        for (core::index2 i = 0; i < y.size(); i++) {
          y[i] = -log(y[i] + pseudocounts_fraction) + log(pseudocounts_fraction);
        }
      }
      logger << utils::LogLevel::FINER << "Registering ff component: " << key << "\n";
      cnt++;
      std::shared_ptr<Interpolate1D<std::vector<double>, core::real, CatmullRomInterpolator<double>>> it = std::make_shared<Interpolate1D<std::vector<double>, core::real, CatmullRomInterpolator<double>>>(x, y, cri);
      std::shared_ptr<Function1D<core::real>> interpolator = it;

      std::shared_ptr<Function1D<core::real>> b_interpolator = std::make_shared<BoundedMFComponent>(interpolator,left_bound,right_bound);
      mf_sp->add_component(key, b_interpolator);
      left_bound = x.front(); // --- reset the bounds to the default values before the next header line is processed
      right_bound = x.back();
    }
  }
  std::vector<std::string> keys = mf_sp->known_distributions();
  if (keys.size() != cnt) {
    logger << utils::LogLevel::SEVERE << ff_file << " file provided " << cnt << "distributions, but distinct loaded only "
        << keys.size() << "\n";
  }

  return mf_sp;
}

const std::vector<std::string> MeanFieldDistributions::known_distributions() const {

  std::vector<std::string> valid_keys;
  std::transform(ff.begin(), ff.end(), std::back_inserter(valid_keys),
      [](const std::map<std::string, EnergyComponent_SP>::value_type& val) {return val.first;});
  return valid_keys;
}

}
}
}

