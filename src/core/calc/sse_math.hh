#ifndef CORE_CALC_sse_math_HH
#define CORE_CALC_sse_math_HH

#include <xmmintrin.h>

#include <core/data/basic/Vec3.hh>
#include <utils/string_utils.hh>

namespace core {
namespace calc {

inline __m128 sse_matrix_times_vector(__m128 c0, __m128 c1, __m128 c2,float a,float b, float c) {

  __m128 va = _mm_set1_ps(a);
  __m128 vc = _mm_mul_ps(va,c0);
  va = _mm_set1_ps(b);
  vc = _mm_add_ps(vc,_mm_mul_ps(va,c1));
  va = _mm_set1_ps(c);
  vc = _mm_add_ps(vc,_mm_mul_ps(va,c2));

  return vc;
}

inline __m128 sse_matrix_times_vector(__m128 c0, __m128 c1, __m128 c2, const core::data::basic::Vec3 & v) {

  return sse_matrix_times_vector(c0,c1,c2,v.x,v.y,v.z);
}

std::ostream & operator<<(std::ostream & out,__m128 r);

}
}

#endif
