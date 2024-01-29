#include <string>
#ifndef UTILS_exit_H
#define UTILS_exit_H

#define EXIT(m) utils::exit_with_error( __FILE__, __LINE__, m )

#define runtime_assert(_Expression, msg) \
        if ( !(_Expression) ) { \
                utils::exit_with_error(__FILE__, __LINE__, msg ); \
        }

namespace utils {

void exit_with_error(const std::string &fname, const int line, const std::string &message);

void exit_OK_with_message(const std::string &message);

}

#endif
