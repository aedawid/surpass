#include <memory>
#include <string>
#include <unordered_map>

#include <core/chemical/Monomer.hh>
#include <core/data/io/Cif.hh>
#include <utils/Logger.hh>
#include <utils/string_utils.hh>

namespace core {
namespace chemical {

using namespace core::data::io;

utils::Logger Monomer::logs = utils::Logger("Monomer");

// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
const Monomer Monomer::ALA(0,'A',"ALA",'P',13,6,0,0.0,0);
const Monomer Monomer::ARG(1,'R',"ARG",'P',27,12,0,0.0,1);
const Monomer Monomer::ASN(2,'N',"ASN",'P',17,9,0,0.0,2);
const Monomer Monomer::ASP(3,'D',"ASP",'P',16,9,0,0.0,3);
const Monomer Monomer::CYS(4,'C',"CYS",'P',14,7,0,0.0,4);
const Monomer Monomer::GLN(5,'Q',"GLN",'P',20,10,0,0.0,5);
const Monomer Monomer::GLU(6,'E',"GLU",'P',19,10,0,0.0,6);
const Monomer Monomer::GLY(7,'G',"GLY",'P',10,5,0,0.0,7);
const Monomer Monomer::HIS(8,'H',"HIS",'P',20,11,0,0.0,8);
const Monomer Monomer::ILE(9,'I',"ILE",'P',22,9,0,0.0,9);
const Monomer Monomer::LEU(10,'L',"LEU",'P',22,9,0,0.0,10);
const Monomer Monomer::LYS(11,'K',"LYS",'P',24,10,0,0.0,11);
const Monomer Monomer::MET(12,'M',"MET",'P',20,9,0,0.0,12);
const Monomer Monomer::PHE(13,'F',"PHE",'P',23,12,0,0.0,13);
const Monomer Monomer::PRO(14,'P',"PRO",'P',17,8,0,0.0,14);
const Monomer Monomer::SER(15,'S',"SER",'P',14,7,0,0.0,15);
const Monomer Monomer::THR(16,'T',"THR",'P',17,8,0,0.0,16);
const Monomer Monomer::TRP(17,'W',"TRP",'P',27,15,0,0.0,17);
const Monomer Monomer::TYR(18,'Y',"TYR",'P',24,13,0,0.0,18);
const Monomer Monomer::VAL(19,'V',"VAL",'P',19,8,0,0.0,19);
const Monomer Monomer::UNK(20,'X',"UNK",'P',5,6,0,0.0,20);
const Monomer Monomer::a(21,'a',"  A",'N',0,6,0,0.0,21);
const Monomer Monomer::c(22,'c',"  C",'N',0,6,0,0.0,22);
const Monomer Monomer::g(23,'g',"  G",'N',0,6,0,0.0,23);
const Monomer Monomer::t(24,'t',"  T",'N',0,6,0,0.0,24);
const Monomer Monomer::u(25,'u',"  U",'N',0,6,0,0.0,25);
const Monomer Monomer::GAP(26,'-',"GAP",'U',0,6,0,0.0,26);
const Monomer Monomer::GPE(27,'-',"GPE",'U',0,6,0,0.0,26);
const Monomer Monomer::UNL(28,'X',"UNL",'U',0,6,0,0.0,28);
const Monomer Monomer::UNG(29,'X',"UNG",'U',0,6,0,0.0,29);

std::unordered_map<std::string,Monomer> Monomer::by_code3 = Monomer::create_map3();
std::unordered_map<char,Monomer> Monomer::by_code1 = Monomer::create_map1();

std::vector<Monomer> & Monomer::by_id() {

  static std::vector<Monomer> by_id( { Monomer::ALA, Monomer::ARG, Monomer::ASN, Monomer::ASP, Monomer::CYS,
      Monomer::GLN, Monomer::GLU, Monomer::GLY, Monomer::HIS, Monomer::ILE, Monomer::LEU, Monomer::LYS, Monomer::MET,
      Monomer::PHE, Monomer::PRO, Monomer::SER, Monomer::THR, Monomer::TRP, Monomer::TYR, Monomer::VAL, Monomer::UNK,
      Monomer::a, Monomer::c, Monomer::g, Monomer::t, Monomer::u, Monomer::GAP, Monomer::GPE, Monomer::UNL, Monomer::UNG });

  return by_id;
}

std::unordered_map<std::string,Monomer> Monomer::create_map3() {

	std::unordered_map<std::string,Monomer> tmp;
	tmp.insert(std::pair<std::string,Monomer>(Monomer::ALA.code3,Monomer::ALA));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::ARG.code3,Monomer::ARG));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::ASN.code3,Monomer::ASN));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::ASP.code3,Monomer::ASP));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::CYS.code3,Monomer::CYS));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::GLN.code3,Monomer::GLN));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::GLU.code3,Monomer::GLU));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::GLY.code3,Monomer::GLY));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::HIS.code3,Monomer::HIS));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::ILE.code3,Monomer::ILE));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::LEU.code3,Monomer::LEU));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::LYS.code3,Monomer::LYS));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::MET.code3,Monomer::MET));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::PHE.code3,Monomer::PHE));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::PRO.code3,Monomer::PRO));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::SER.code3,Monomer::SER));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::THR.code3,Monomer::THR));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::TRP.code3,Monomer::TRP));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::TYR.code3,Monomer::TYR));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::VAL.code3,Monomer::VAL));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::UNK.code3,Monomer::UNK));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::a.code3,Monomer::a));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::c.code3,Monomer::c));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::g.code3,Monomer::g));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::t.code3,Monomer::t));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::u.code3,Monomer::u));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::GAP.code3,Monomer::GAP));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::GPE.code3,Monomer::GPE));
	tmp.insert(std::pair<std::string,Monomer>(Monomer::UNL.code3,Monomer::UNL));

	return tmp;
}

std::unordered_map<char,Monomer>Monomer:: create_map1() {

	std::unordered_map<char,Monomer> tmp;
	tmp.insert(std::pair<char,Monomer>(Monomer::ALA.code1,Monomer::ALA));
	tmp.insert(std::pair<char,Monomer>(Monomer::ARG.code1,Monomer::ARG));
	tmp.insert(std::pair<char,Monomer>(Monomer::ASN.code1,Monomer::ASN));
	tmp.insert(std::pair<char,Monomer>(Monomer::ASP.code1,Monomer::ASP));
	tmp.insert(std::pair<char,Monomer>(Monomer::CYS.code1,Monomer::CYS));
	tmp.insert(std::pair<char,Monomer>(Monomer::GLN.code1,Monomer::GLN));
	tmp.insert(std::pair<char,Monomer>(Monomer::GLU.code1,Monomer::GLU));
	tmp.insert(std::pair<char,Monomer>(Monomer::GLY.code1,Monomer::GLY));
	tmp.insert(std::pair<char,Monomer>(Monomer::HIS.code1,Monomer::HIS));
	tmp.insert(std::pair<char,Monomer>(Monomer::ILE.code1,Monomer::ILE));
	tmp.insert(std::pair<char,Monomer>(Monomer::LEU.code1,Monomer::LEU));
	tmp.insert(std::pair<char,Monomer>(Monomer::LYS.code1,Monomer::LYS));
	tmp.insert(std::pair<char,Monomer>(Monomer::MET.code1,Monomer::MET));
	tmp.insert(std::pair<char,Monomer>(Monomer::PHE.code1,Monomer::PHE));
	tmp.insert(std::pair<char,Monomer>(Monomer::PRO.code1,Monomer::PRO));
	tmp.insert(std::pair<char,Monomer>(Monomer::SER.code1,Monomer::SER));
	tmp.insert(std::pair<char,Monomer>(Monomer::THR.code1,Monomer::THR));
	tmp.insert(std::pair<char,Monomer>(Monomer::TRP.code1,Monomer::TRP));
	tmp.insert(std::pair<char,Monomer>(Monomer::TYR.code1,Monomer::TYR));
	tmp.insert(std::pair<char,Monomer>(Monomer::VAL.code1,Monomer::VAL));
	tmp.insert(std::pair<char,Monomer>(Monomer::UNK.code1,Monomer::UNK));
	tmp.insert(std::pair<char,Monomer>(Monomer::a.code1,Monomer::a));
	tmp.insert(std::pair<char,Monomer>(Monomer::c.code1,Monomer::c));
	tmp.insert(std::pair<char,Monomer>(Monomer::g.code1,Monomer::g));
	tmp.insert(std::pair<char,Monomer>(Monomer::t.code1,Monomer::t));
	tmp.insert(std::pair<char,Monomer>(Monomer::u.code1,Monomer::u));
	tmp.insert(std::pair<char,Monomer>(Monomer::GAP.code1,Monomer::GAP));
	tmp.insert(std::pair<char,Monomer>(Monomer::GPE.code1,Monomer::GPE));
	tmp.insert(std::pair<char,Monomer>(Monomer::UNL.code1,Monomer::UNL));

	return tmp;
}

//  1020 20 ? U 255 0  2.0 AQS
Monomer::Monomer(const std::string & line) {

  const char * str = line.c_str();
  id = utils::to_int(str);
  parent_id = utils::to_int(str + 5);
  code1 = line[9];
  type = line[11];
  n_atoms = utils::to_int(str + 13);
  n_heavy_atoms = n_atoms;
//  n_heavy_atoms = utils::to_int(str + 13);
  is_ambiguous = (line[17]=='0') ? false : true;
  charge = utils::to_double(str + 19);
  code3 = std::string(line.substr(24, 3));
  if (type == 'N') code1 = tolower(code1);
}


Monomer::Monomer(const core::index2 id,const char code1,const std::string code3,const char type,const unsigned char n_atoms,const unsigned char n_heavy_atoms,
    const bool ambig_flag,const core::real charge,const core::index2 parent) :
    id(id), parent_id(parent), code1(code1), type(type), n_atoms(n_atoms), n_heavy_atoms(n_heavy_atoms), is_ambiguous(
        ambig_flag), charge(charge), code3(code3) {
}

Monomer::Monomer(const core::index2 id,const char code1,const char* code3,const char type,const unsigned char n_atoms,const unsigned char n_heavy_atoms,
    const bool ambig_flag,const core::real charge,const core::index2 parent) :
    id(id), parent_id(parent), code1(code1), type(type), n_atoms(n_atoms), n_heavy_atoms(n_heavy_atoms), is_ambiguous(
        ambig_flag), charge(charge), code3(code3) {
}

Monomer::Monomer(const Monomer &m) :
    Monomer(m.id, m.code1, m.code3, m.type, m.n_atoms, m.n_heavy_atoms, m.is_ambiguous, m.charge, m.parent_id) {
}


std::ostream& operator<< (std::ostream &out, const Monomer &m) {

	out << utils::string_format("%4d %c %3s %c %d %4.1f %4d",m.id, m.code1, m.code3.c_str(), m.type, m.is_ambiguous, m.charge, m.parent_id);
    return out;
}

std::ostream& operator<< (std::ostream &out, const Monomer *m) {

	out << utils::string_format("%4d %c %3s %c %d %4.1f %4d",m->id, m->code1, m->code3.c_str(), m->type, m->is_ambiguous, m->charge, m->parent_id);
    return out;
}

}
}
