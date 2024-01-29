#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <mutex>

#include <utils/Logger.hh>
#include <utils/LogManager.hh>
#include <utils/string_utils.hh>

namespace utils {

std::mutex Logger::mtx;
std::thread::id Logger::blocking_thread;

bool Logger::is_logable(LogLevel level) const { return (level <= LogManager::log_level); }

Logger &operator <<(Logger &logger, const char c) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << c;

  return logger;
}


Logger &operator <<(Logger &logger, const float value) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << value;

  return logger;
}

Logger &operator <<(Logger &logger, const char* message) {

  if (logger.recent_level <= LogManager::log_level) {
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << message;
  }

  return logger;
}

Logger &operator <<(Logger &logger, const size_t number) {

  if (logger.recent_level <= LogManager::log_level) {
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << number;
  }

  return logger;
}


Logger &operator <<(Logger &logger, const unsigned long long number) {

  if (logger.recent_level <= LogManager::log_level) {
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << number;
  }

  return logger;
}

Logger &operator <<(Logger &logger, const int number) {

  if (logger.recent_level <= LogManager::log_level) {
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << number;
  }

  return logger;
}

Logger &operator <<(Logger &logger, const double value) {

  if (logger.recent_level <= LogManager::log_level) {
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << value;
  }

  return logger;
}

Logger &operator <<(Logger &logger, const std::string & message) {

  if (logger.recent_level <= LogManager::log_level) {
    if (logger.blocking_thread != std::this_thread::get_id()) {
      std::lock_guard<std::mutex> lock(Logger::mtx);
      if (!LogManager::is_muted(logger.module_name_)) logger.sink << message;
    } else if (!LogManager::is_muted(logger.module_name_)) logger.sink << message;
  }

  return logger;
}

Logger &operator <<(Logger &logger, const LogLevel level) {

  if (LogManager::is_muted(logger.module_name_)) return logger;

  if (logger.blocking_thread != std::this_thread::get_id()) {
    std::lock_guard<std::mutex> lock(Logger::mtx);

    logger.recent_level = level;
    if (logger.recent_level <= LogManager::log_level) {
      if (LogManager::use_colors())
        logger.sink << log_level_names_colored.at(level) << TEXT_BOLD << logger.module_name_ << TEXT_RESET << " ";
      else
        logger.sink << log_level_names.at(level) << logger.module_name_ << " ";
    }
  } else {
    logger.recent_level = level;
    if (logger.recent_level <= LogManager::log_level) {
      if(LogManager::use_colors())
        logger.sink << log_level_names_colored.at(level) << TEXT_BOLD << logger.module_name_ << TEXT_RESET << " ";
      else
        logger.sink << log_level_names.at(level) <<  logger.module_name_ << " ";
    }
  }
  return logger;
}

}
