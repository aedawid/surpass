#include <cstdio>
#include <cstdlib>

#ifndef _WIN32
#include <unistd.h>
#endif

#include <iostream>

#include <utils/exit.hh>

namespace utils {

void exit_with_error(const std::string &fname, const int line, const std::string &message) {

  if (isatty(fileno(stdout))) std::cerr << "\x1b[0m\x1b[1m\x1b[31m";
  if (!message.empty()) std::cerr << std::endl << "ERROR: " << message << std::endl;
  std::cerr << "ERROR:: Exit from: " << fname << " line: " << line << std::endl;
  if (isatty(fileno(stdout))) std::cerr << "\x1b[0m";
  std::cerr.flush();
  std::exit(1);
}

void exit_OK_with_message(const std::string &message) {

  if (isatty(fileno(stdout))) std::cerr << "\x1b[0m\x1b[1m\x1b[31m";
  if (!message.empty()) std::cerr << std::endl << message << std::endl;
  if (isatty(fileno(stdout))) std::cerr << "\x1b[0m";
  std::cerr.flush();
  std::exit(0);
}

}
