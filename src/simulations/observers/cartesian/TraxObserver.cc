
#include <iostream>
#include <fstream>

#include <core/data/basic/Vec3.hh>
#include <core/data/basic/Vec3Cubic.hh>

#include <utils/string_utils.hh>

#include <simulations/observers/cartesian/TraxObserver.hh>

namespace simulations {
namespace observers {
namespace cartesian {

template<typename C>
TraxObserver<C>::TraxObserver(const systems::ResidueChain<C> &observed_object, const core::index1 max_atom_in_residue,
             const std::string & out_fname) : ObserveEvaluators(std::make_shared<std::stringstream >()),
             observed_object(observed_object),n_frames(0), atoms_per_line(max_atom_in_residue), out_fname(out_fname) {

  outstream = std::make_shared<std::ofstream>(out_fname);
  evaluator_stream = std::static_pointer_cast<std::stringstream>(ObserveEvaluators::output_stream());
}

template<typename C>
bool TraxObserver<C>::observe() {

  ++cnt;
  if(!ObserverInterface::trigger->operator()()) return false;

  ++n_frames;
  core::data::basic::Vec3 cm;

  for (core::index4 i = 0; i < observed_object.n_atoms; i++) cm += observed_object.coordinates[i];
  cm /= double(observed_object.n_atoms);
  std::string tag = utils::string_format("%02d%05d\n",observed_object.system_id(),cnt);
  (*outstream) << utils::string_format("# %9.3f %9.3f %9.3f",cm.x,cm.y,cm.z);
  if (energy_ != nullptr) (*outstream) << *energy_;
  ObserveEvaluators::observe();
  std::string s (evaluator_stream->str());
  s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());
  (*outstream) << s << " " << tag;
  evaluator_stream->str("");
  for(core::index2 ires=0;ires < observed_object.count_residues();++ires) {
    auto atoms_range = observed_object.atoms_for_residue(ires);
    char restype = core::chemical::Monomer::get(core::index2(observed_object.residue_type(atoms_range.first_atom))).code1;
    (*outstream) << utils::string_format("%3d %c ",ires+1,restype);

    for (auto iatom = atoms_range.first_atom; iatom <= atoms_range.last_atom; ++iatom)
      (*outstream) << utils::string_format("%7.3f %7.3f %7.3f ", observed_object[iatom].x - cm.x,
        observed_object[iatom].y - cm.y, observed_object[iatom].z - cm.z);

    for (auto iatom = 0; iatom < atoms_per_line - atoms_range.last_atom + atoms_range.first_atom - 1; ++iatom)
      (*outstream) << utils::string_format("%7.3f %7.3f %7.3f ", 0,0,0);
    (*outstream) << tag;
  }
  outstream->flush();

  return true;
}

template<typename C>
std::string TraxObserver<C>::header_string() const {
  std::stringstream s;
  s << " center-x  center-y  center-z  ";
  if (energy_ != nullptr) s << energy_->header_string();
  s << ObserveEvaluators::header_string() << " tag";

  return s.str();
}

template class TraxObserver<core::data::basic::Vec3>;
template class TraxObserver<core::data::basic::Vec3Cubic>;

}
}
}
