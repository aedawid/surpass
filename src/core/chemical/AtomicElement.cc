#include <core/chemical/AtomicElement.hh>
#include <utils/Logger.hh>

namespace core {
namespace chemical {

utils::Logger AtomicElement::logger("AtomicElement");

const AtomicElement AtomicElement::DUMMY = AtomicElement(0,0.0,"Du","Dummy");
const AtomicElement AtomicElement::HYDROGEN = AtomicElement(1,1.008,"H","Hydrogen");
const AtomicElement AtomicElement::DEUTERIUM = AtomicElement(1,2.013553,"D","Deuterium");
const AtomicElement AtomicElement::HELIUM = AtomicElement(2,4.002602,"He","Helium");
const AtomicElement AtomicElement::LITHIUM = AtomicElement(3,6.941000,"Li","Lithium");
const AtomicElement AtomicElement::BERYLLIUM = AtomicElement(4,9.012182,"Be","Beryllium");
const AtomicElement AtomicElement::BORON = AtomicElement(5,10.811000,"B","Boron");
const AtomicElement AtomicElement::CARBON = AtomicElement(6,12.010700,"C","Carbon");
const AtomicElement AtomicElement::NITROGEN = AtomicElement(7,14.006700,"N","Nitrogen");
const AtomicElement AtomicElement::OXYGEN = AtomicElement(8,15.999400,"O","Oxygen");
const AtomicElement AtomicElement::FLUORINE = AtomicElement(9,18.998403,"F","Fluorine");
const AtomicElement AtomicElement::NEON = AtomicElement(10,20.179700,"Ne","Neon");
const AtomicElement AtomicElement::SODIUM = AtomicElement(11,22.989769,"Na","Sodium");
const AtomicElement AtomicElement::MAGNESIUM = AtomicElement(12,24.305000,"Mg","Magnesium");
const AtomicElement AtomicElement::ALUMINIUM = AtomicElement(13,26.981539,"Al","Aluminium");
const AtomicElement AtomicElement::SILICON = AtomicElement(14,28.085500,"Si","Silicon");
const AtomicElement AtomicElement::PHOSPHORUS = AtomicElement(15,30.973762,"P","Phosphorus");
const AtomicElement AtomicElement::SULFUR = AtomicElement(16,32.065000,"S","Sulfur");
const AtomicElement AtomicElement::CHLORINE = AtomicElement(17,35.453000,"Cl","Chlorine");
const AtomicElement AtomicElement::POTASSIUM = AtomicElement(19,39.098300,"K","Potassium");
const AtomicElement AtomicElement::ARGON = AtomicElement(18,39.948000,"Ar","Argon");
const AtomicElement AtomicElement::CALCIUM = AtomicElement(20,40.078000,"Ca","Calcium");
const AtomicElement AtomicElement::SCANDIUM = AtomicElement(21,44.955912,"Sc","Scandium");
const AtomicElement AtomicElement::TITANIUM = AtomicElement(22,47.867000,"Ti","Titanium");
const AtomicElement AtomicElement::VANADIUM = AtomicElement(23,50.941500,"V","Vanadium");
const AtomicElement AtomicElement::CHROMIUM = AtomicElement(24,51.996100,"Cr","Chromium");
const AtomicElement AtomicElement::MANGANESE = AtomicElement(25,54.938045,"Mn","Manganese");
const AtomicElement AtomicElement::IRON = AtomicElement(26,55.845000,"Fe","Iron");
const AtomicElement AtomicElement::NICKEL = AtomicElement(28,58.693400,"Ni","Nickel");
const AtomicElement AtomicElement::COBALT = AtomicElement(27,58.933195,"Co","Cobalt");
const AtomicElement AtomicElement::COPPER = AtomicElement(29,63.546000,"Cu","Copper");
const AtomicElement AtomicElement::ZINC = AtomicElement(30,65.409000,"Zn","Zinc");
const AtomicElement AtomicElement::GALLIUM = AtomicElement(31,69.723000,"Ga","Gallium");
const AtomicElement AtomicElement::GERMANIUM = AtomicElement(32,72.640000,"Ge","Germanium");
const AtomicElement AtomicElement::ARSENIC = AtomicElement(33,74.921600,"As","Arsenic");
const AtomicElement AtomicElement::SELENIUM = AtomicElement(34,78.960000,"Se","Selenium");
const AtomicElement AtomicElement::BROMINE = AtomicElement(35,79.904000,"Br","Bromine");
const AtomicElement AtomicElement::KRYPTON = AtomicElement(36,83.798000,"Kr","Krypton");
const AtomicElement AtomicElement::RUBIDIUM = AtomicElement(37,85.467800,"Rb","Rubidium");
const AtomicElement AtomicElement::STRONTIUM = AtomicElement(38,87.620000,"Sr","Strontium");
const AtomicElement AtomicElement::YTTRIUM = AtomicElement(39,88.905850,"Y","Yttrium");
const AtomicElement AtomicElement::ZIRCONIUM = AtomicElement(40,91.224000,"Zr","Zirconium");
const AtomicElement AtomicElement::NIOBIUM = AtomicElement(41,92.906000,"Nb","Niobium");
const AtomicElement AtomicElement::MOLYBDENUM = AtomicElement(42,95.940000,"Mo","Molybdenum");
const AtomicElement AtomicElement::TECHNETIUM = AtomicElement(43,0.000000,"Tc","Technetium");
const AtomicElement AtomicElement::RUTHENIUM = AtomicElement(44,101.070000,"Ru","Ruthenium");
const AtomicElement AtomicElement::RHODIUM = AtomicElement(45,102.905000,"Rh","Rhodium");
const AtomicElement AtomicElement::PALLADIUM = AtomicElement(46,106.420000,"Pd","Palladium");
const AtomicElement AtomicElement::SILVER = AtomicElement(47,107.868200,"Ag","Silver");
const AtomicElement AtomicElement::CADMIUM = AtomicElement(48,112.411000,"Cd","Cadmium");
const AtomicElement AtomicElement::INDIUM = AtomicElement(49,114.818000,"In","Indium");
const AtomicElement AtomicElement::TIN = AtomicElement(50,118.710000,"Sn","Tin");
const AtomicElement AtomicElement::ANTIMONY = AtomicElement(51,121.760000,"Sb","Antimony");
const AtomicElement AtomicElement::TELLURIUM = AtomicElement(52,127.600000,"Te","Tellurium");
const AtomicElement AtomicElement::IODINE = AtomicElement(53,126.904000,"I","Iodine");
const AtomicElement AtomicElement::XENON = AtomicElement(54,131.293000,"Xe","Xenon");
const AtomicElement AtomicElement::CAESIUM = AtomicElement(55,132.905452,"Cs","Caesium");
const AtomicElement AtomicElement::BARIUM = AtomicElement(56,137.327000,"Ba","Barium");
const AtomicElement AtomicElement::LANTHANUM = AtomicElement(57,138.905470,"La","Lanthanum");
const AtomicElement AtomicElement::CERIUM = AtomicElement(58,140.116000,"Ce","Cerium");
const AtomicElement AtomicElement::PRASEODYMIUM = AtomicElement(59,140.907650,"Pr","Praseodymium");
const AtomicElement AtomicElement::NEODYMIUM = AtomicElement(60,144.242000,"Nd","Neodymium");
const AtomicElement AtomicElement::PROMETHIUM = AtomicElement(61,0.000000,"Pm","Promethium");
const AtomicElement AtomicElement::SAMARIUM = AtomicElement(62,150.360000,"Sm","Samarium");
const AtomicElement AtomicElement::EUROPIUM = AtomicElement(63,151.964000,"Eu","Europium");
const AtomicElement AtomicElement::GADOLINIUM = AtomicElement(64,157.250000,"Gd","Gadolinium");
const AtomicElement AtomicElement::TERBIUM = AtomicElement(65,158.925350,"Tb","Terbium");
const AtomicElement AtomicElement::DYSPROSIUM = AtomicElement(66,162.500000,"Dy","Dysprosium");
const AtomicElement AtomicElement::HOLMIUM = AtomicElement(67,164.930000,"Ho","Holmium");
const AtomicElement AtomicElement::ERBIUM = AtomicElement(68,167.259000,"Er","Erbium");
const AtomicElement AtomicElement::THULIUM = AtomicElement(69,168.934210,"Tm","Thulium");
const AtomicElement AtomicElement::YTTERBIUM = AtomicElement(70,173.040000,"Yb","Ytterbium");
const AtomicElement AtomicElement::LUTETIUM = AtomicElement(71,174.967000,"Lu","Lutetium");
const AtomicElement AtomicElement::HAFNIUM = AtomicElement(72,178.490000,"Hf","Hafnium");
const AtomicElement AtomicElement::TANTALUM = AtomicElement(73,180.947880,"Ta","Tantalum");
const AtomicElement AtomicElement::TUNGSTEN = AtomicElement(74,183.840000,"W","Tungsten");
const AtomicElement AtomicElement::RHENIUM = AtomicElement(75,186.207000,"Re","Rhenium");
const AtomicElement AtomicElement::OSMIUM = AtomicElement(76,190.230000,"Os","Osmium");
const AtomicElement AtomicElement::IRIDIUM = AtomicElement(77,192.217000,"Ir","Iridium");
const AtomicElement AtomicElement::PLATINUM = AtomicElement(78,195.084000,"Pt","Platinum");
const AtomicElement AtomicElement::GOLD = AtomicElement(79,196.966569,"Au","Gold");
const AtomicElement AtomicElement::MERCURY = AtomicElement(80,200.590000,"Hg","Mercury");
const AtomicElement AtomicElement::THALLIUM = AtomicElement(81,204.383300,"Tl","Thallium");
const AtomicElement AtomicElement::LEAD = AtomicElement(82,207.200000,"Pb","Lead");
//const AtomicElement AtomicElement::URANIUM = AtomicElement(92,238.028910,"U","Uranium");
const AtomicElement AtomicElement::BISMUTH = AtomicElement(83,208.980400,"Bi","Bismuth");

const AtomicElement AtomicElement::periodic_table[] = {DUMMY, HYDROGEN, HELIUM, LITHIUM, BERYLLIUM, BORON, CARBON, NITROGEN, OXYGEN, FLUORINE, NEON, SODIUM, MAGNESIUM, ALUMINIUM, SILICON, PHOSPHORUS, SULFUR, CHLORINE, ARGON, POTASSIUM, CALCIUM, SCANDIUM, TITANIUM, VANADIUM, CHROMIUM, MANGANESE, IRON, COBALT, NICKEL, COPPER, ZINC, GALLIUM, GERMANIUM, ARSENIC, SELENIUM, BROMINE, KRYPTON, RUBIDIUM, STRONTIUM, YTTRIUM, ZIRCONIUM, NIOBIUM, MOLYBDENUM, TECHNETIUM, RUTHENIUM, RHODIUM, PALLADIUM, SILVER, CADMIUM, INDIUM, TIN, ANTIMONY, TELLURIUM, IODINE, XENON, CAESIUM, BARIUM, LANTHANUM, CERIUM, PRASEODYMIUM, NEODYMIUM, PROMETHIUM, SAMARIUM, EUROPIUM, GADOLINIUM, TERBIUM, DYSPROSIUM, HOLMIUM, ERBIUM, THULIUM, YTTERBIUM, LUTETIUM, HAFNIUM, TANTALUM, TUNGSTEN, RHENIUM, OSMIUM, IRIDIUM, PLATINUM, GOLD, MERCURY, THALLIUM, LEAD, BISMUTH};
const std::unordered_map<std::string,AtomicElement> AtomicElement::elements_by_symbol = AtomicElement::create_map();

const AtomicElement & AtomicElement::by_mass(const core::real mass) {

  auto hit = std::lower_bound(periodic_table, (periodic_table) + elements_by_symbol.size(), mass,
    [](const AtomicElement &e1, const core::real m) { return e1.mass < m; });
  core::real d1 = ((*hit).mass - mass);
  d1 = d1*d1;
  core::real d2 = ((*(hit-1)).mass - mass);
  d2 *= d2;
  return (d1<d2) ? (*hit) : *(hit-1);
}

const AtomicElement & AtomicElement::by_symbol(const std::string &symbol) {

  try {
    return elements_by_symbol.at(symbol);
  } catch (std::out_of_range &e) {
    if (symbol.size() == 2) {
      std::string s2(symbol);
      s2[1] = std::tolower(s2[1]);
      try {
        return elements_by_symbol.at(s2);
      } catch (std::out_of_range &e) {
        logger << utils::LogLevel::WARNING << "can't find atomic element for the string:>" << s2 << "<\n";
        return DUMMY;
      }
    }
    logger << utils::LogLevel::WARNING << "can't find atomic element for the string:>" << symbol << "<\n";
    return DUMMY;
  }
}

std::ostream& operator<< (std::ostream &out, const AtomicElement &e) {

    out << e.z<<" "<<e.symbol << " : " << e.name;
    return out;
}

}
}

