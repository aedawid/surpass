#include <memory>
#include <string>
#include <vector>
#include <cstdlib> // for std::exit
#include <iostream>
#include <stdexcept>

#include <core/index.hh>
#include <core/SURPASSenvironment.hh>
#include <core/SURPASSversion.hh>

#include <utils/options/Option.hh>
#include <utils/options/OptionParser.hh>
#include <utils/options/input_options.hh>
#include <utils/string_utils.hh>

namespace utils {
namespace options {

utils::Logger OptionParser::logger = utils::Logger("OptionParser");

const std::string & OptionParser::value_string(const Option & o) const { return value_string(o.name); }

const std::string & OptionParser::value_string(const std::string & o) const {

  if (option_called.size() == 0) {
    logger << utils::LogLevel::SEVERE << "\nNo command-line options called... Use -h to see what's available\n";
    return empty;
  }
  const auto & which_option = registered_keys.find(o);
  if (which_option == registered_keys.end()) {
    throw std::runtime_error("\nAttempt to read " + o + " option which was not used at command line");
  }
  return option_values[(which_option->second)];
}

bool OptionParser::was_used(const Option & o) const {

  if (option_called.size() == 0) {
    logger << utils::LogLevel::SEVERE << "\nNo command-line options called... Use -h to see what's available\n";
    return 0;
  }
  const auto & which_option = registered_keys.find(o.name);
  if (which_option == registered_keys.end()) {
    logger << utils::LogLevel::SEVERE
        << "Testing option " + o.name + " which was not declared by the program as a valid option!\n";
    return false;
  }

  return option_called[which_option->second];
}

bool OptionParser::was_used(const std::string &o) const {

  if (option_called.size() == 0) {
    logger << utils::LogLevel::SEVERE << "\nNo command-line options called... Use -h to see what's available\n";
    return 0;
  }
  const auto & which_option = registered_keys.find(o);
  if (which_option == registered_keys.end()) {
    logger << utils::LogLevel::SEVERE
           << "Testing option " + o + " which was not declared by the program as a valid option!\n";
    return false;
  }

  return option_called[which_option->second];
}

bool OptionParser::parse_cmdline(const core::index2 argc, const char *args[],bool exit_if_empty) {

  logger << utils::LogLevel::INFO << "\n" << core::SURPASSversion().to_string();

  // ---------- If user has not given any option, print help and exit
  if (exit_if_empty) if (argc < 2) print_help(std::cerr);

  program_name_ = args[0];

  bool db_path_from_environment = false;
  // ---------- Otherwise, parse the options
  for (core::index2 i = 0; i < argc; ++i) {
    if (args[i][0] != '-') continue;
    const std::string o(args[i]);
    std::vector<std::string> tokens;
    if (o.find('=') == std::string::npos) tokens.push_back(o);
    else tokens = utils::split(o, '=');
    if ((tokens[0] == std::string("-verbose")) || (tokens[0] == std::string("-v"))) {
      if( tokens.size() != 2 ) 
        std::cerr << "Option -verbose (or -v) requires a parameter, e.g. -verbose=WARNING to see only program warnings\n";
      else LogManager::get().set_level(tokens[1]);
    }
    if ((tokens[0] == std::string("-in:database")) || (tokens[0] == std::string("-d"))) {
      core::SURPASSenvironment::surpass_db_path(tokens[1]);
    } else {
      std::string path = core::get_env_var("SURPASS_DATA_DIR");
      if(path.size()>1) {
        db_path_from_environment = true;
        core::SURPASSenvironment::surpass_db_path(path);
        logger << utils::LogLevel::INFO << "SURPASS DB path extracted from a shell variable: "
            << core::SURPASSenvironment::surpass_db_path() << "\n";
      }
    }

    if (tokens[0] == std::string("-mute")) {
      for (const std::string & l : utils::split(tokens[1], ','))
        utils::LogManager::mute(l);
    }

    if ((tokens[0] == std::string("-help")) || (tokens[0] == std::string("-h")))  print_help(std::cerr);

    if ((tokens[0] == std::string("-help:markdown")) || (tokens[0] == std::string("-help::markdown")) || (tokens[0] == std::string("-hmd")))  print_help_markdown(std::cerr);

    if (tokens[0] == std::string("-examples")) {
      print_examples(std::cerr);
    }
  }
  option_called.resize(registered_options.size(), 0);
  option_values.resize(registered_options.size(),"");
  for (core::index2 i = 0; i < registered_options.size(); ++i)
    option_values[i] = registered_options[i].default_value_string;


  for (core::index2 i = 0; i < argc; ++i) {
    if (args[i][0] != '-') continue;
    std::string o(args[i]);
    std::vector<std::string> tokens;
    if (o.find('=') == std::string::npos) tokens.push_back(o);
    else tokens = utils::split(o, '=');
    utils::replace_substring(tokens[0], "::", ":");
    utils::replace_substring(tokens[0], "--", "-");

    auto which_option = registered_keys.find(tokens[0]);
    logger << utils::LogLevel::FINE << "Option " << tokens[0] << " has been seen in the command line\n";

    // ---------- Print log-message about all keys recognized by this program
    if (which_option == registered_keys.end()) {
      logger << utils::LogLevel::SEVERE << "\nRequested option " << tokens[0]
          << " was not declared by the program as a valid option!\n";
      std::for_each(registered_keys.begin(),registered_keys.end(),[](const std::pair<std::string,size_t> & p){ std::cerr << p.first << "\n"; });
      throw std::invalid_argument(
          "Requested option " + tokens[0] + " was not declared by the program as a valid option!\n");
    }

    size_t option_index = which_option->second;
    Option & op = registered_options[option_index];

    option_called[option_index] = true;
    if (tokens.size() == 2) {
      if(tokens[1].size()==0) {
        std::string msg = utils::string_format("Empty string given as the value for an option %s\n",tokens[0].c_str());
        logger << utils::LogLevel::SEVERE << msg;
        throw std::invalid_argument(msg);
      }
      option_values[option_index] = tokens[1];

      logger << utils::LogLevel::FINE << "Value " << option_values[option_index] << " associated with option "
          << op.get_name() << "\n";
    } else {
      if (tokens.size() > 2)
        logger << utils::LogLevel::SEVERE << "Can't parse give value (too many '=') : " << tokens[1] << "\n";
      logger << utils::LogLevel::SEVERE << "Can't parse give value (too many '=') : " << tokens[1] << "\n";

      if(value_mandatory[option_index]) {
        std::string msg = utils::string_format("Option %s must provided a value (e.g. %s=some_value)\n",tokens[0].c_str(),tokens[0].c_str());
        logger << utils::LogLevel::SEVERE << msg;
        throw std::invalid_argument(msg);
      }
    }
  }

  // ---------- check if all mandatory options have shown up
  bool was_ok = true;
  for (Option &op : registered_options) {
    if (op.mandatory()) {
      logger << utils::LogLevel::SEVERE << "The mandatory option: " << op.get_name() << " was not used\n";
      was_ok = false;
    }
  }

  if(is_registered(db_path)&&db_path_from_environment) {
    const size_t k = registered_keys[db_path.name];
    option_values[k] = core::SURPASSenvironment::surpass_db_path();
    option_called[k]=true;
  }

  already_parsed = true;

  if (!was_ok) {
//		print_help(std::cerr);
    return false;
  }

  return true;
}

void OptionParser::inject_option(Option & o, const std::string & value) {

  unsigned short n_opt = registered_keys.size();
  if (!is_registered(o)) {
    register_option(o);
    for (unsigned short i = n_opt; i < registered_keys.size(); ++i) {
      option_called.push_back(true);
      option_values.push_back(value);
    }
  } else {
    auto which_option = registered_keys.find(o.name);
    option_values[which_option->second] = value;
    option_called[which_option->second] = true;
  }

  logger << utils::LogLevel::INFO << "Option " << o.name << " injected with value '" << value << "'\n";
}

void OptionParser::print_help(std::ostream & where) {

  where << "\n" << program_info_ << "\n";

  core::index2 longest = 0;
  for (const Option &op : registered_options)
    if (op.get_name().length() > longest) longest = op.get_name().length();
  where << "\n";
  std::string line_padding(longest+2,' ');
  std::string fmt = utils::string_format("%%%ds :", longest);
  for (const Option &op : registered_options) {
    std::vector<std::string> words = utils::split(op.get_info(),' ');
    std::string info = utils::format_paragraph(words,"",line_padding,100,100-longest-3);
    where << utils::string_format(fmt, op.get_name().c_str())<<info<<"\n";
  }
  if (groups_of_examples.size() > 0)
    where << TEXT_BLUE << "\n\nUse:\n\n\t" << program_name_ << " -examples\n\nto see examples !\n\n" << TEXT_RESET;
  else where << "\n";

  std::exit(0);
}

void OptionParser::print_help_markdown(std::ostream & where) {

  where << "\n" << program_info_ << "\n";

  core::index2 longest = 0;
  for (const Option &op : registered_options)
    if (op.get_name().length() > longest) longest = op.get_name().length();
  where << "\n";
  std::string line_padding(longest+2,' ');

  where << "| ---:|:---|\n";
  std::string fmt = utils::string_format("|%%%ds | %%s|\n", longest);
  for (const Option &op : registered_options) {
    where << utils::string_format(fmt, op.get_name().c_str(),op.get_info().c_str());
  }
  where << "\nHere are a few examples:\n\n";

  for (size_t i = 0; i < groups_of_examples.size(); i++) {
    where  << example_group_names[i] << "\n\n";
    for (core::index2 j = 0; j <  groups_of_examples[i].size(); ++j) {
      where << "+ " << groups_of_examples[i][j].first << "\n\n```\n" << groups_of_examples[i][j].second << "\n```\n\n";
    }
    where << "\n\n";
  }

  std::exit(0);
}

void OptionParser::print_examples(std::ostream & where,core::index2 screen_width) {

  where << "\n";
  for (size_t i = 0; i < groups_of_examples.size(); i++) {
    core::index2 n = (screen_width - example_group_names[i].size()) / 2;
    where << TEXT_BLUE<< std::string(n - 1, '-') << " " << example_group_names[i] << " " << std::string(n - 1, '-')
                             << TEXT_RESET<< "\n";
    for (core::index2 j = 0; j <  groups_of_examples[i].size(); ++j) {
      where << (j + 1) << ") " << groups_of_examples[i][j].first << "\n" << groups_of_examples[i][j].second << "\n\n";
    }
    where << "\n";
  }
  std::exit(0);
}

bool OptionParser::is_registered(const Option & o) const {

  return is_registered(o.name);
}

bool OptionParser::is_registered(const std::string & full_name) const {

  if (registered_keys.find(full_name) == registered_keys.end()) return false;
  return true;
}

std::string OptionParser::build_cmdline() {

  std::stringstream out;
  out << program_name_<<" ";
  for(auto o :registered_keys) {
    if(option_called[o.second]) {
      out << o.first;
      if (option_values[o.second] != "") out << "=" << option_values[o.second] << " ";
      else out << " ";
    }
  }
  return out.str();
}

bool OptionParser::register_option(Option &o) {

  if (is_registered(o)) return false;
  registered_options.push_back(o);
  size_t index = registered_options.size() - 1;
  std::string long_name = o.name;
  utils::replace_substring(long_name, "::", ":");
  utils::replace_substring(long_name, "--", "-");
  registered_keys[o.name] = index;
  registered_keys[long_name] = index;
  registered_keys[o.short_name] = index;
  value_mandatory.push_back(false);
  logger << LogLevel::FINER << "Registering " << o.name << " " << o.short_name << "\n";

  return true;
}


void OptionParser::add_example(core::index2 which_group,const std::string& example_info,const std::string & example_cmd) {

  if(already_parsed) throw std::logic_error("Command line has been already parsed! No more examples may be added!");

  if (which_group >= example_group_names.size()) throw std::out_of_range("Examples' group index too high!");

  groups_of_examples[which_group].push_back(std::pair<std::string,std::string>(example_info,example_cmd));
}

void OptionParser::add_example(core::index2 which_group,const char* example_info,const std::string & example_cmd) {

  if(already_parsed) throw std::logic_error("Command line has been already parsed! No more examples may be added!");

  add_example(which_group, std::string(example_info), std::string(example_cmd));
}

void VerbosityLevel::execute() { utils::LogManager::set_level(utils::options::OptionParser::get().value_string("-v")); }

void PrintHelp::execute() { utils::options::OptionParser::get().print_help(std::cerr); }

void PrintMarkdownHelp::execute() { utils::options::OptionParser::get().print_help_markdown(std::cerr); }

}
}

