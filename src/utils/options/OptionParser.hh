/** @file OptionParser.hh
 *
 * Provides OptionParser object and a few staple option objects.
 * More option object are defined in other header files:
 *	    - utils/options/select_options.hh
 *	    - utils/options/Option.hh
 *	    - utils/options/align_options.hh
 *	    - utils/options/calc_options.hh
 *	    - utils/options/input_options.hh
 *	    - utils/options/output_options.hh
 */
#ifndef UTILS_OPTIONS_OptionParser_HH
#define UTILS_OPTIONS_OptionParser_HH

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <exception>

#include <core/index.hh>
#include <utils/options/Option.hh>
#include <utils/Logger.hh>
#include <utils/string_utils.hh>

namespace utils {
namespace options {

/** @brief Command line parser for the programs.
 *
 * The parser takes <code>argc</code> and <code>argv</code> arguments that <code>main()</code> function received
 * from operating system and parses their content.
 *
 * This is a small example that shows how to use the option system:
 * @include ex_OptionParser.cc
 */
class OptionParser {
public:

  /** @brief Returns the singleton instance of a command line parser.
   *
   * @param program_name - name of the program
   * @return singleton reference to the parser
   */
  static OptionParser &get(const std::string &program_name = "__UNINITIALIZED_") {
    static OptionParser manager(program_name);
    return manager;
  }


/** @name Register one or more options for the parser.
 *
 * Only registered options will be parsed at runtime. The use of unregistered option will stop the program.
 */
///@{
  bool register_option(Option &o);

  bool register_option(Option &o1, Option &o2) { return register_option(o1) && register_option(o2); }

  bool register_option(Option &o1, Option &o2, Option &o3) { return register_option(o1) && register_option(o2, o3); }

  bool register_option(Option &o1, Option &o2, Option &o3, Option &o4) {
    return register_option(o1, o2) && register_option(o3, o4);
  }

  bool register_option(Option &o1, Option &o2, Option &o3, Option &o4, Option &o5) {
    return register_option(o1, o2) && register_option(o3, o4) && register_option(o5);
  }

  bool register_option(Option &o1, Option &o2, Option &o3, Option &o4, Option &o5, Option &o6) {
    return register_option(o1, o2) && register_option(o3, o4) && register_option(o5, o6);
  }

  bool register_option(Option &o1, Option &o2, Option &o3, Option &o4, Option &o5, Option &o6, Option &o7) {
    return register_option(o1, o2) && register_option(o3, o4) && register_option(o5, o6) && register_option(o7);
  }

  bool register_option(Option &o1, Option &o2, Option &o3, Option &o4, Option &o5, Option &o6, Option &o7, Option &o8) {
    return register_option(o1, o2) && register_option(o3, o4) && register_option(o5, o6) && register_option(o7, o8);
  }
///@}

  /** @brief returns true if a given option has been already registered for the parser
   *
   * @param o - option object
   * @return true when already registered, false otherwise
   */
  bool is_registered(const Option &o) const;

  /** @brief returns true if a given option has been already registered for the parser
   *
   * @param full_name - full name of the option
   * @return true when already registered, false otherwise
   */
  bool is_registered(const std::string &full_name) const;

  /** @brief Tells a program a certain option was used from command line, even if it wasn't
   *
   * @param o - an option
   * @param value - value assigned to this option
   */
  void inject_option(Option &o, const std::string &value);

  /** @brief Parses the command line.
   *
   * The key method of this class should be called at the beginning of the program but already after all the relevant
   * options have been already registered. It parses the command line, extracts values and assigns them to their options.
   *
   * @param argc - starndard parameter of the <code>main()</code> function
   * @param args  - starndard parameter of the <code>main()</code> function
   * @param exit_if_empty - should the program quit if the command line was empty?
   * @return true when parsed successfully
   */
  bool parse_cmdline(const core::index2 argc, const char *args[], bool exit_if_empty = true);

  /// Sets the help message for the program
  void program_info(const std::string &msg) { program_info_ = msg; }

  /** @brief Writes a help message to a stream.
   *
   * @param where - a stream the text is sent to
   */
  void print_help(std::ostream &where);

  /** @brief Writes a help message in <strong>markdown</strong> format to a stream.
   *
   * @param where - a stream the text is sent to
   */
  void print_help_markdown(std::ostream &where);

  /** @brief Returns the value (as a string) assigned to the given option object
   *
   * @param o - option object
   * @return value from command line - as a string. If no value was given at runtime, then the default value is returned.
   * If this is also not defined, empty string is returned
   */
  const std::string &value_string(const Option &o) const;

  /** @brief Returns the value (as a string) assigned to the given option object
   *
   * @param o - option name
   * @return same as <code>std::string &value_string(const Option &)</code>
   */
  const std::string &value_string(const std::string &o) const;

  /// Returns true if an option was used at a command line
  bool was_used(const Option &o) const;

  /// Returns true if an option was used at a command line
  bool was_used(const std::string &o) const;

  /** @brief Returns the name of the program that's running.
   *
   * The name is actually <code>argv[0]</code> element extracted from the cmdline that was parsed
   */
  const std::string &program_name() const { return program_name_; }

  /** @brief Re-creates the command line string from the options user used and their respective values
   *
   * @return a command line string that reproduces the behavior of a program.
   * The order of options in the returned sting may be different from the original call.
   */
  std::string build_cmdline();

  /** @brief Define that an option must be acompanied with a value
   *
   * @param o - an option string
   * @param must_have_value - true or false
   */
  void must_have_value(const Option &o, const bool must_have_value) {
    value_mandatory[registered_keys.find(o.name)->second] = must_have_value;
  }

  /** @brief Define that an option must be acompanied with a value
   *
   * @param o - an option string
   */
  bool must_have_value(const Option &o) { return value_mandatory[registered_keys.find(o.name)->second]; }

  /** @brief Creates a new group of examples.
   *
   * The group is empty i.e. there is no single example in it. Examples should be added with
   * <code> add_example()</code> method
   *
   * @param group_name - name of this group of examples - it will be printed on the screen when program's help is invoked
   * @return index of the new example group; use this index to add examples to the group
   */
  inline core::index2 add_examples_group(const std::string &group_name) {

    example_group_names.push_back(group_name);
    groups_of_examples.push_back(std::vector<std::pair<std::string, std::string>>());
    return example_group_names.size() - 1;
  }

  /** @brief Add a new example to a group.
   *
   * @param which_group - index of the respective group of examples
   * @param example_info - string that describes this example, e.g. what does it do - just a single line of text
   * @param example_cmd - example command
   */
  void add_example(core::index2 which_group, const std::string &example_info, const std::string &example_cmd);

  /** @brief Add a new example to a group.
   *
   * @param which_group - index of the respective group of examples
   * @param example_info - string that describes this example, e.g. what does it do - just a single line of text
   * @param example_cmd - example command
   */
  void add_example(core::index2 which_group, const char *example_info, const std::string &example_cmd);

  /** @brief Prints all examples on the screen
   *
   * @param where - output stream (typically std::cout)
   * @param screen_width - with of the text column
   */
  void print_examples(std::ostream &where, core::index2 screen_width = 80);

private:
  static utils::Logger logger;
  std::vector<Option> registered_options;
  bool already_parsed;
  std::string program_name_;
  bool is_initialized = false;

  std::vector<std::vector<std::pair<std::string, std::string>>> groups_of_examples;
  std::vector<std::string> example_group_names;

  std::unordered_map<std::string, size_t> registered_keys;
  std::vector<std::string> option_values;
  std::vector<bool> option_called;
  std::vector<bool> value_mandatory;
  const std::string empty = "";
  std::string program_info_;

  OptionParser(const std::string &program_name) {
    program_name_ = program_name;
    already_parsed = false;
  }

  OptionParser(OptionParser const &);

  void operator=(OptionParser const &);
};

/** @brief Returns a value associated with a given option.
 *
 * When no value was specified at the command line, the default value is returned
 *
 * @param o - option object
 * @tparam T  - cast the value to this type
 * @return option value
 */
template<typename T>
inline T option_value(const Option &o) {
  return utils::from_string<T>(OptionParser::get().value_string(o));
}

/** @brief Returns a value associated with a given option.
 *
 * When no value was specified at the command line, the default value is returned
 *
 * @param o - option name
 * @tparam T  - cast the value to this type
 * @return option value
 */
template<typename T>
inline T option_value(const std::string &o) {
  return utils::from_string<T>(OptionParser::get().value_string(o));
}

/** @brief Returns values associated with a given option.
 *
 * When no value was specified at the command line, the default value string is parsed and its content returned
 *
 * @param o - option object
 * @tparam T  - cast the value to this type
 * @tparam C  - container type for the values, e.g. <code>std::vector<T> </code>
 * @return option value
 */
template<typename C, typename T>
inline C option_value(const Option &o) {

  C container;
  std::vector<std::string> tokens = utils::split(OptionParser::get().value_string(o), ',');
  for (const std::string & t : tokens)
    container.push_back(utils::from_string<T>(t));

  return container;
}

/** @brief Returns values associated with a given option.
 *
 * When no value was specified at the command line, the default value string is parsed and its content returned
 *
 * @param o - option name
 * @tparam T  - cast the value to this type
 * @tparam C  - container type for the values, e.g. <code>std::vector<T> </code>
 * @return option value
 */
template<typename C, typename T>
inline C option_value(const std::string &o) {

  C container;
  std::vector<std::string> tokens = utils::split(OptionParser::get().value_string(o), ',');
  for (const std::string & t : tokens)
    container.push_back(utils::from_string<T>(t));

  return container;
}

/** @brief Returns a value associated with a given option.
 *
 * This methods overrides the default value option with a given value.
 *
 * @param o - option object
 * @tparam T  - cast the value to this type
 * @return  value associated with a given option. If no value was specified at the command line, the given
 * <code>default_value</code> value is casted to the type <code>T</code> and returned
 */
template<typename T>
inline T option_value(const Option &o, const T default_value) {

  if ((o.was_used()) && (OptionParser::get().value_string(o) != ""))
    return utils::from_string<T>(OptionParser::get().value_string(o));

  else return default_value;
}

/** @brief Returns a value associated with a given option.
 *
 * This methods overrides the default value option with a given value.
 *
 * @param o - option name
 * @tparam T  - cast the value to this type
 * @return  value associated with a given option. If no value was specified at the command line, the given
 * <code>default_value</code> value is casted to the type <code>T</code> and returned
 */
template<typename T>
inline T option_value(const std::string &o, const T default_value) {

  if ((OptionParser::get().was_used(o)) && (OptionParser::get().value_string(o) != ""))
    return utils::from_string<T>(OptionParser::get().value_string(o));

  else return default_value;
}

/// Explicit speciation for std::string type
template<>
inline std::string option_value<std::string>(const Option &o) {
  return OptionParser::get().value_string(o);
}

/// Explicit speciation for std::string type
template<>
inline std::string option_value<std::string>(const std::string &o) {
  return OptionParser::get().value_string(o);
}

/** @brief Stores values associated with a given option in a given vector.
 *
 * The vector is not cleared before copying, so the values are appended to it.
 * When no value was specified at the command line, the default value string is parsed and its content returned
 *
 * @param o - option name
 * @param destination - vector where the values associated with this option will be appended
 * @tparam T  - cast the value to this type
 * @return option value
 */
template<typename T>
std::vector<T> &  option_value(const Option &o, std::vector<T> &destination) {

  if (o.was_used()) split(option_value<std::string>(o), destination, ',');

  return destination;
}

/// Ask for silencing one or several loggers (to limit output on stderr)
static Option mute("-mute", "-mute", "mute a particular logger; more than one logger names may be specified when separated with a comma");

/// Ask to execute program in several threads
static Option n_threads("-t", "-n_threads", "asks an executable to run calculations in more than one thread", "1");

/// Change verbosity level of a program
class VerbosityLevel : public Option {
public:
  VerbosityLevel() : Option("-v", "-verbose", "set the verbosity level", "", false) {}
  virtual void execute();
}static verbose;


/// Print examples for a program
static Option show_examples("-examples", "-examples", "print just the examples, not the full help message", "", false);

/// Print program's help message in markdown format
class PrintMarkdownHelp : public Option {
public:
  PrintMarkdownHelp() : Option("-hmd", "-help:markdown", "print help message in markdown format (to be included in the on-line documentation)", "", false) {}
  virtual void execute();
}static markdown_help;

/// Print program's help message
class PrintHelp : public Option {
public:
  PrintHelp() : Option("-h", "-help", "print help message", "", false) {}
  virtual void execute();
}static help;
}
}
/**
 * @example ex_OptionParser.cc
 */
#endif

