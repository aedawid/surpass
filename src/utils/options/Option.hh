#ifndef UTILS_OPTIONS_Option_HH
#define UTILS_OPTIONS_Option_HH

#include <iostream>
#include <stdexcept>
#include <string>
#include <memory>

namespace utils {
namespace options {

class OptionParser;

/// Object that corresponds to a single option (flag) that may be specified at command line for any program
class Option {
public:

  const std::string short_name;   ///< a shortcut for this option, e.g. <tt>-h</tt> for <tt>-help</tt>
  const std::string name;         ///< option name i.e. the text given at command line, e.g.  <tt>-help</tt> for <tt>-in::pdb</tt>
  const std::string info;         ///< help message for this option
  const std::string default_value_string; ///< default value for this option - as a string (will be parsed to its final type later)

  /** Creates a new option object
   *
   * @param short_name - short (single-character) name for this option, e.g. <code>-h</code> or <code>-n</code>
   * @param name - long name for this option e.g. <code>-help</code>
   * @param info - description string for this option that will be printed in a help message
   * @param default_value_string - default value for this option
   * @param if_mandatory - if true, the option is mandatory.
   */
  Option(const std::string &short_name, const std::string &name, const std::string &info,
    const std::string &default_value_string = "", const bool if_mandatory = false) :
    short_name(short_name),name(name), info(info), default_value_string(default_value_string),mandatory_(if_mandatory) {}

  /// Virtual destructor
  virtual ~Option() { }

  /// Returns the long name of this option
  inline const std::string & get_name() const { return name; }

  /// Returns the description string assiciated with this option
  inline const std::string & get_info() const { return info; }

  /// Returns the short name of this option
  inline const std::string & get_short_name() const { return short_name; }

  /// True if this option was used at the command line
  bool was_used() const;

  virtual void execute() {}

  /// Declare this option as mandatory
  void mandatory(const bool if_mandatory) { mandatory_ = if_mandatory; }

  /// Returns true if this option was declared mandatory
  bool mandatory() { return mandatory_; }

private:
  bool mandatory_;        ///< if set to true, a program will require this option to show up at command line
  std::shared_ptr<OptionParser> manager;
};

}
}
#endif
