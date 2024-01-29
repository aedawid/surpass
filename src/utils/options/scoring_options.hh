/** \file scoring_options.hh
 * @brief defines option objects that control energy evaluations
 */
#ifndef UTILS_OPTIONS_scoring_options_HH
#define UTILS_OPTIONS_scoring_options_HH

#include <utils/options/Option.hh>

namespace utils {
namespace options {

/** @name Options that control which energy function will be used for scoring
 */
///@{
static Option scfx_config("-x", "-in:scfx", "provide a config file that defines energy function");
static Option cabs_bb_go("-cabs_bb_go", "-scfx:cabs_bb_go", "use default CABS-bb with additional Go-like term; requires reference structure.");
static Option cabs_bb("-cabs_bb", "-scfx:cabs_bb", "use default CABS-bb (CABS with explicit backbone) energy for scoring");
static Option cabs("-cabs", "-scfx:cabs_bb", "use default CABS energy for scoring");
static Option surpass("-surpass", "-scfx:surpass", "use default SURPASS energy for scoring");
///@}

}
}

#endif
