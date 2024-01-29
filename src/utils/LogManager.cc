#include <utils/Logger.hh>
#include <utils/LogManager.hh>

namespace utils {

bool LogManager::use_colors_ = true;
bool LogManager::use_colors_setup_ = false;
LogLevel LogManager::log_level = LogLevel::INFO;
utils::Logger LogManager::logger = utils::Logger("LogLevel");
std::unordered_set<std::string> LogManager::muted;

bool LogManager::use_colors() {

  if (!use_colors_setup_) {
    if (!isatty(fileno(stderr))) use_colors(false);
    else use_colors(true);
    use_colors_setup_ = true;
  }
  return use_colors_;
}

}
