#include <core/data/basic/Vec3.hh>

#include <core/calc/structural/transformations/Rototranslation.hh>
#include <core/data/structural/PdbAtom.hh>

using core::data::basic::Vec3;

namespace core {
namespace calc {
namespace structural {
namespace transformations {

Rototranslation_SP local_coordinates_three_atoms(const Vec3 &a1, const Vec3 &a2, const Vec3 &a3) {

  Rototranslation_SP t = std::make_shared<Rototranslation>();

  Vec3 nToCa(a2);// a1 -> a2 vector
  nToCa -= a1;
  nToCa.norm();

  Vec3 caToC(a3);// a2 -> a3 vector
  caToC -= a2;
  caToC.norm();

  Vec3 tz(nToCa);
  tz -= caToC;
  tz.norm();

  Vec3 tx(nToCa);
  tx += caToC;
  tx.norm();
  Vec3 ty;
  data::basic::cross_product(nToCa, caToC, ty);
  ty.norm();

  t->tr_before(a2);
  t->tr_after(0.0, 0.0, 0.0);
  t->rot_x(tx.x, tx.y, tx.z);
  t->rot_y(ty.x, ty.y, ty.z);
  t->rot_z(tz.x, tz.y, tz.z);

  return t;
}

Rototranslation_SP local_coordinates_three_atoms(const core::data::structural::Residue &r,
                                                 const std::vector<std::string> &atom_names) {

  return local_coordinates_three_atoms(*r.find_atom_safe(atom_names[0]), *r.find_atom_safe(atom_names[1]),
    *r.find_atom_safe(atom_names[2]));
}

void local_coordinates_three_atoms(const Vec3 &a1, const Vec3 &a2, const Vec3 &a3, Rototranslation & r) {

  Vec3 nToCa(a2);// a1 -> a2 vector
  nToCa -= a1;
  nToCa.norm();

  Vec3 caToC(a3);// a2 -> a3 vector
  caToC -= a2;
  caToC.norm();


  Vec3 tz(nToCa);
  tz -= caToC;
  tz.norm();

  Vec3 tx(nToCa);
  tx += caToC;
  tx.norm();
  Vec3 ty;
  data::basic::cross_product(nToCa, caToC, ty);
  ty.norm();

  r.tr_before_.set(a2);
  r.tr_after_.set(0.0, 0.0, 0.0);
  r.rot_x_.set(tx.x, tx.y, tx.z);
  r.rot_y_.set(ty.x, ty.y, ty.z);
  r.rot_z_.set(tz.x, tz.y, tz.z);
}


Vec3 euler_angles(const Vec3 &row_x, const Vec3 &row_y, const Vec3 &row_z) {
  Vec3 out(-asin(row_z.x), atan2(row_z.y, row_z.z), atan2(row_y.x, row_x.x));
  return out;
}

Vec3 euler_angles(const Rototranslation &r) {
  return euler_angles(r.rot_x(), r.rot_y(), r.rot_z());
}

Vec3 euler_angles(const Rototranslation &r1, const Rototranslation &r2) {

  Vec3 v1, v2, v3;

  v1.x = r1.rot_x().x * r2.rot_x().x + r1.rot_x().y * r2.rot_y().x + r1.rot_x().z * r2.rot_z().x;
  v1.y = r1.rot_x().x * r2.rot_x().y + r1.rot_x().y * r2.rot_y().y + r1.rot_x().z * r2.rot_z().y;
  v1.z = r1.rot_x().x * r2.rot_x().z + r1.rot_x().y * r2.rot_y().z + r1.rot_x().z * r2.rot_z().z;

  v2.x = r1.rot_y().x * r2.rot_x().x + r1.rot_y().y * r2.rot_y().x + r1.rot_y().z * r2.rot_z().x;
  v2.y = r1.rot_y().x * r2.rot_x().y + r1.rot_y().y * r2.rot_y().y + r1.rot_y().z * r2.rot_z().y;
  v2.z = r1.rot_y().x * r2.rot_x().z + r1.rot_y().y * r2.rot_y().z + r1.rot_y().z * r2.rot_z().z;

  v3.x = r1.rot_z().x * r2.rot_x().x + r1.rot_z().y * r2.rot_y().x + r1.rot_z().z * r2.rot_z().x;
  v3.y = r1.rot_z().x * r2.rot_x().y + r1.rot_z().y * r2.rot_y().y + r1.rot_z().z * r2.rot_z().y;
  v3.z = r1.rot_z().x * r2.rot_x().z + r1.rot_z().y * r2.rot_y().z + r1.rot_z().z * r2.rot_z().z;

  return euler_angles(v1, v2, v3);
}

}
}
}
}
