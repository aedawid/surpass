#include <string>
#include <fstream>
#include <stdexcept>

#include <core/SURPASSenvironment.hh>
#include <utils/io_utils.hh>
#include <utils/options/OptionParser.hh>
#include <utils/options/input_options.hh>

namespace core {

std::string SURPASSenvironment::surpass_db_path_ = "";
bool SURPASSenvironment::db_path_env_tested_ = false;
utils::Logger SURPASSenvironment::logger("SURPASSenvironment");

const std::string get_env_var(const std::string & key) {
  char * val;
  val = getenv( key.c_str() );
  std::string retval = "";
  if (val != NULL) {
    retval = val;
  }

  return retval;
}
    
void SURPASSenvironment::set_db_path_from_environment() {
  
  std::string path = core::get_env_var("SURPASS_DATA_DIR");
  SURPASSenvironment::db_path_env_tested_ = true;
  if(path.size()>1) {
    SURPASSenvironment::surpass_db_path(path);
    logger << utils::LogLevel::INFO << "SURPASS DB path extracted from a shell variable: "
        << core::SURPASSenvironment::surpass_db_path() << "\n";
    }
  }

const std::string SURPASSenvironment::from_file_or_db(const std::string & fname) {

  using utils::options::option_value;
  bool if_found;
  std::string fname2;
  if(!db_path_env_tested_)
    set_db_path_from_environment();
  std::string the_path_now = surpass_db_path_;
  if (utils::options::db_path.was_used())
    the_path_now = option_value<std::string>(utils::options::db_path);
  if (the_path_now.size()>0) {
    if_found = utils::find_file(fname, fname2, the_path_now);
    if (if_found) logger << utils::LogLevel::FILE << "file " << fname2 << " found in SURPASS DB location: "
        << the_path_now << "\n";
    else logger << utils::LogLevel::FILE << "Can't find file " << fname2 << " in SURPASS DB location: "
        << option_value<std::string>(utils::options::db_path) << "\n";
  } else {
    logger << utils::LogLevel::FILE << "DB path not set. Use a proper command line option or set the SURPASS_DATA_DIR environmenal (i.e. shell) variable.\nFile not found: " << fname<<"\n";
    if_found = utils::find_file(fname, fname2, ".");
  }
  if (!if_found) throw std::runtime_error("Can't locate " + fname);
  return fname2;
}

const std::string SURPASSenvironment::from_file_or_db(const std::string & fname, const std::string & db_location) {

  using utils::options::option_value;
  bool if_found;
  std::string fname2;
  if(!db_path_env_tested_)
    set_db_path_from_environment();
  std::string the_path_now = surpass_db_path_;
  if (utils::options::db_path.was_used())
    the_path_now = option_value<std::string>(utils::options::db_path);

  if (the_path_now.size()>0) {
    const std::string db_path = utils::join_paths(the_path_now, db_location);
    if_found = utils::find_file(fname, fname2, db_path);
    if (if_found) logger << utils::LogLevel::FILE << "file " << fname2 << " found in SURPASS DB location: " << db_path
        << "\n";
    else logger << utils::LogLevel::FILE << "Can't find file " << fname2 << " in SURPASS DB location: " << db_path
        << "\n";
  } else {
      logger << utils::LogLevel::FILE << "DB path not set. Use a proper command line option or set the SURPASS_DATA_DIR environmenal (i.e. shell) variable.\nFile not found: " << fname<<"\n";
    if_found = utils::find_file(fname, fname2, ".");
  }
  if (!if_found)
    throw std::runtime_error("Can't locate " + fname);

  return fname2;
}

void SURPASSenvironment::surpass_db_path(const std::string & new_path) {
  surpass_db_path_ = new_path;
  logger << utils::LogLevel::FILE << "Setting the SURPASS database path to " << surpass_db_path_ << "\n";
}

const std::string & SURPASSenvironment::surpass_db_path() {

  if((surpass_db_path_=="")&&(!db_path_env_tested_))
    set_db_path_from_environment();

  return surpass_db_path_;
}

}

