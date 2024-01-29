#ifndef SIMULATIONS_FORCEFIELDS_ForceFieldConfig_HH
#define SIMULATIONS_FORCEFIELDS_ForceFieldConfig_HH

#include <unordered_map>
#include <vector>
#include <string>
#include <utils/Logger.hh>


namespace simulations {
namespace forcefields {

class ForceFieldConfig {
public:

  static std::vector<std::string> known_substitutions;

  ForceFieldConfig(const std::string & config_as_txt) : logger("ForceFieldConfig"), cfg(config_as_txt) {}

  void input_pdb(const std::string & input_pdb_fname) { set("${INPUT_PDB}",input_pdb_fname); }

  const std::string & input_pdb() const { return substitutions.at("${INPUT_PDB}"); }

  void input_ss2(const std::string & input_ss2_fname) { set("${INPUT_SS2}", input_ss2_fname); }

  const std::string & input_ss2() const { return substitutions.at("${INPUT_SS2}"); }

  void native_pdb(const std::string & native_pdb_fname) { set("${NATIVE_PDB}", native_pdb_fname); }

  const std::string & native_pdb() const { return substitutions.at("${NATIVE_PDB}"); }

  const std::string &  substitute();

  const std::string &  str() const { return cfg; }

private:
  utils::Logger logger;
  std::string cfg;
  std::unordered_map<std::string,std::string> substitutions;

  void set(const std::string & key, const std::string & value);
};

} // ~ simulations
} // ~ forcefields
#endif
