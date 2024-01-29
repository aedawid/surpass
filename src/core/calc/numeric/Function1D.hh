#ifndef CORE_CALC_NUMERIC_Function1D_HH
#define CORE_CALC_NUMERIC_Function1D_HH

#include <memory>

#include <core/real.hh>

namespace core {
namespace calc {
namespace numeric {

/** @brief Pure virtual definition of a 1D function.
 *
 * This interface defines a one-argument call operator.
 */
template <typename E>
class Function1D {
public:

  /** @brief Call a 1D function
   * @param x - function's argument
   */
  virtual E operator()(const E x) const = 0;
};

}
}
}

#endif
