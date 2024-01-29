#ifndef UTILS_Logger_HH
#define UTILS_Logger_HH

#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>
#include <mutex>
#include <thread>

#include <utils/string_utils.hh>

namespace utils {

enum class LogLevel {
	CRITICAL = 0,
	SEVERE = 1,
	WARNING = 2,
	FILE = 3,
	INFO = 4,
	FINE = 5,
	FINER = 6,
	FINEST = 7,
};

static const std::pair<LogLevel, std::string> log_level_pairs[] = { //
		std::pair<LogLevel, std::string>(LogLevel::CRITICAL, "[CRITICAL] "), //
		std::pair<LogLevel, std::string>(LogLevel::SEVERE,   "  [SEVERE] "), //
		std::pair<LogLevel, std::string>(LogLevel::WARNING,  " [WARNING] "), //
		std::pair<LogLevel, std::string>(LogLevel::INFO,     "    [INFO] "), //
		std::pair<LogLevel, std::string>(LogLevel::FINE,     "    [FINE] "), //
		std::pair<LogLevel, std::string>(LogLevel::FINER,    "   [FINER] "), //
		std::pair<LogLevel, std::string>(LogLevel::FINEST,   "  [FINEST] "), //
		std::pair<LogLevel, std::string>(LogLevel::FILE,     "    [FILE] ") //
		};

static const std::map<LogLevel, std::string> log_level_names(log_level_pairs,
    log_level_pairs + sizeof(log_level_pairs) / sizeof(log_level_pairs[0]));

static const std::pair<LogLevel, std::string> log_level_pairs_colored[] = { //
    std::pair<LogLevel, std::string>(LogLevel::CRITICAL, std::string(TEXT_RED)+std::string("[CRITICAL] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::SEVERE,   std::string(TEXT_RED)+std::string("  [SEVERE] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::WARNING,  std::string(TEXT_YELLOW)+std::string(" [WARNING] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::INFO,     std::string(TEXT_YELLOW)+std::string("    [INFO] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::FINE,     std::string(TEXT_WHITE)+std::string("    [FINE] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::FINER,    std::string(TEXT_WHITE)+std::string("   [FINER] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::FINEST,   std::string(TEXT_WHITE)+std::string("  [FINEST] ")+std::string(TEXT_RESET)), //
    std::pair<LogLevel, std::string>(LogLevel::FILE,     std::string(TEXT_WHITE)+std::string("    [FILE] ")+std::string(TEXT_RESET)) //
    };

static const std::map<LogLevel, std::string> log_level_names_colored(log_level_pairs_colored,
    log_level_pairs_colored + sizeof(log_level_pairs_colored) / sizeof(log_level_pairs_colored[0]));


class Logger {
private:
  static std::mutex mtx;
  static std::thread::id blocking_thread;

public:

	Logger(const std::string & module_name) :
			module_name_(module_name), recent_level(LogLevel::INFO) {
	}

	inline const std::string & module_name() const { return module_name_; }

	virtual ~Logger() {}

	bool is_logable(LogLevel level) const;

	friend Logger &operator <<(Logger &logger, const LogLevel level);

	friend Logger &operator <<(Logger &logger, const char* message);

	friend Logger &operator <<(Logger &logger, const std::string & message);

	friend Logger &operator <<(Logger &logger, const double value);

	friend Logger &operator <<(Logger &logger, const int number);

	friend Logger &operator <<(Logger &logger, const size_t number);

	friend Logger &operator <<(Logger &logger, const unsigned long long number);

  friend Logger &operator <<(Logger &logger, const float value);

  friend Logger &operator <<(Logger &logger, const char c);

  template <class ...Args>
	void log(const Args& ...args) {

    std::lock_guard<std::mutex> lock(mtx);
    blocking_thread = std::this_thread::get_id();
    log_one(args...);
    sink.flush();
	}

private:
	const std::string  module_name_;
	std::ostream &sink = std::cerr;
	LogLevel recent_level;

	template <class A0, class ...Args>
  void log_one(const A0& a0, const Args& ...args) {
	  *this  << a0;
	  log_one(args...);
  }

	inline void log_one() {}
};

}
#endif /* UTILS_LOGGER_HH */
