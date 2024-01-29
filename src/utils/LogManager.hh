#ifndef UTILS_LogManager_HH
#define UTILS_LogManager_HH

#include <unistd.h>

#include <unordered_set>
#include <mutex>

#include <utils/Logger.hh>

namespace utils {

class LogManager {
public:
  static LogManager & get() {
    static LogManager manager;
    return manager;
  }

  /// Returns true if loggers print in color
  static bool use_colors();

  /// Switch ON/OFF color output
  static void use_colors(bool value) { use_colors_ = value;}

  static void mute(const std::string & logger_name) {
    logger <<LogLevel::INFO << "muted channel "<<logger_name<<"\n";
    muted.insert(logger_name);
  }

  static bool is_muted(const std::string & logger_name) {
    return (muted.find(logger_name) != muted.end());
  }

  static void set_level(LogLevel & new_level) {
    manager.log_level = new_level;
  }

  static void set_level(const std::string & new_level_name) {

    std::map<LogLevel, std::string>::const_iterator it;
    std::string n(new_level_name);
    for (auto pos = n.begin(); pos != n.end(); ++pos) *pos = toupper(*pos);
    for (it = log_level_names.begin(); it != log_level_names.end(); it++) {
      if(it->second.find(n)!=it->second.npos) {
//      if (it->second.compare(new_level_name) == 0) {
        manager.log_level = it->first;
        logger << it->first << "Logging level set to " << new_level_name << "\n";
        return;
      }
    }
    logger << LogLevel::WARNING << "Unknown log-level: " << new_level_name << ", the value not set. Known levels are:";
    for (it = log_level_names.begin(); it != log_level_names.end(); it++)
      logger << " " << it->second;
    logger << "\n";
  }

  static void CRITICAL() {
    manager.log_level = LogLevel::CRITICAL;
  }
  static void SEVERE() {
    manager.log_level = LogLevel::SEVERE;
  }
  static void WARNING() {
    manager.log_level = LogLevel::WARNING;
  }
  static void FILE() {
    manager.log_level = LogLevel::FILE;
  }
  static void INFO() {
    manager.log_level = LogLevel::INFO;
  }
  static void FINE() {
    manager.log_level = LogLevel::FINE;
  }
  static void FINER() {
    manager.log_level = LogLevel::FINER;
  }
  static void FINEST() {
    manager.log_level = LogLevel::FINEST;
  }

private:
  static bool use_colors_;
  static bool use_colors_setup_;
  static const LogManager manager;
  std::ostream &sink;
  static LogLevel log_level;
  static utils::Logger logger;
  static std::unordered_set<std::string> muted;

  LogManager() : sink(std::cerr) {

    if (!isatty(fileno(stderr))) use_colors(false);
    else use_colors(true);
  }

  LogManager(LogManager const&);
  void operator=(LogManager const&);

  friend Logger;

  friend Logger &operator <<(Logger &logger, const LogLevel level);

  friend Logger &operator <<(Logger &logger, const char* message);

  friend Logger &operator <<(Logger &logger, const std::string & message);

  friend Logger &operator <<(Logger &logger, const double value);

  friend Logger &operator <<(Logger &logger, const int number);

  friend Logger &operator <<(Logger &logger, const size_t number);

  friend Logger &operator <<(Logger &logger, const unsigned long long number);

  friend Logger &operator <<(Logger &logger, const float value);

  friend Logger &operator <<(Logger &logger, const char c);
};

}

#endif // ~ UTILS_LogManager_HH
