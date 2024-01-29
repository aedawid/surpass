#include <unordered_map>
#include <vector>
#include <string>

#include <utils/string_utils.hh>

#include <simulations/forcefields/ForceFieldConfig.hh>

namespace simulations {
namespace forcefields {

std::vector<std::string> ForceFieldConfig::known_substitutions = {"${INPUT_PDB}", "${INPUT_SS2}", "${NATIVE_PDB}"};

const std::string & ForceFieldConfig::substitute() {

  for(auto sub : substitutions) {
    utils::replace_substring(cfg, sub.first, sub.second);
    logger << utils::LogLevel::INFO << "Substituting " << sub.first << " with " << sub.second << "\n";
  }

  return cfg;
}

void ForceFieldConfig::set(const std::string & key, const std::string & value) {

  if(value=="")
    logger << utils::LogLevel::SEVERE << "Substitution string is empty for the keyword: " << key << "\n";

  if (std::find(known_substitutions.cbegin(), known_substitutions.cend(), key) == known_substitutions.cend())
    logger << utils::LogLevel::SEVERE << "Unknown substitution keyword: " << key << "\n";

  substitutions[key] = value;
  substitute();
}
} // ~ simulations
} // ~ forcefields
