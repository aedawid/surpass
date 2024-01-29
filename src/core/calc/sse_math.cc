#include <xmmintrin.h>

#include <core/calc/sse_math.hh>

namespace core {
namespace calc {

std::ostream & operator<<(std::ostream & out,__m128 r) {

  out << utils::string_format("%f %f %f %f",r[0],r[1],r[2],r[3]);

  return out;
}

}
}
