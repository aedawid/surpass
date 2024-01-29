#ifndef CORE_SURPASSenvironment_HH
#define CORE_SURPASSenvironment_HH

#include <string>
#include <fstream>

#include <utils/Logger.hh>

namespace core {

/** @brief Returns a shell variable value
 *
 * @param key - the variable name, e.g. SURPASS_DATA_DIR
 * @return the variable's value
 */
const std::string get_env_var(const std::string & key);

/** @brief defines global constants and environment for SURPASS's programs
 */
class SURPASSenvironment {
public:

  /**@brief Returns the directory where the SURPASS keeps his stuff
   *
   * The path should be defined by the -in:database command line option
   * @return the path to SURPASS's parameters
   */
  static const std::string & surpass_db_path();

  /**@brief Sets the new value of the directory name where the SURPASS keeps his stuff
   *@param new_path - the new location
   */
  static void surpass_db_path(const std::string & new_path);

  static void set_db_path_from_environment();
    
  /**@brief tries to open a file from SURPASS database; local directory is also tested.
   *
   * @param fname - name of the file to be opened
   * @return file name with a properly completed path
   */
  static const std::string from_file_or_db(const std::string & fname);

  /** @brief Looks for a file in the database; local directory is also tested.
   *
   * @param fname - name of the file to be opened
   * @param db_location - relative path within the SURPASS database where the file should be located
   * @return file name with a properly completed path
   */
  static const std::string from_file_or_db(const std::string & fname, const std::string & db_location);

private:
  static utils::Logger logger;
  /// path pointing to the main directory with a database
  static std::string surpass_db_path_;
  static bool db_path_env_tested_;
};

}

#endif
