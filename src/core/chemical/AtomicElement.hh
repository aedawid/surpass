#ifndef CORE_CHEMICAL_AtomicElement_HH
#define CORE_CHEMICAL_AtomicElement_HH

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <stdexcept>

#include <core/index.hh>
#include <core/real.hh>
#include <utils/Logger.hh>

namespace core {
namespace chemical {

using core::real;

class AtomicElement {

public:
  static const AtomicElement periodic_table[];
  static const AtomicElement DUMMY;
  static const AtomicElement HYDROGEN;
  static const AtomicElement DEUTERIUM;
  static const AtomicElement HELIUM;
  static const AtomicElement LITHIUM;
  static const AtomicElement BERYLLIUM;
  static const AtomicElement BORON;
  static const AtomicElement CARBON;
  static const AtomicElement NITROGEN;
  static const AtomicElement OXYGEN;
  static const AtomicElement FLUORINE;
  static const AtomicElement NEON;
  static const AtomicElement SODIUM;
  static const AtomicElement MAGNESIUM;
  static const AtomicElement ALUMINIUM;
  static const AtomicElement SILICON;
  static const AtomicElement PHOSPHORUS;
  static const AtomicElement SULFUR;
  static const AtomicElement CHLORINE;
  static const AtomicElement POTASSIUM;
  static const AtomicElement ARGON;
  static const AtomicElement CALCIUM;
  static const AtomicElement SCANDIUM;
  static const AtomicElement TITANIUM;
  static const AtomicElement VANADIUM;
  static const AtomicElement CHROMIUM;
  static const AtomicElement MANGANESE;
  static const AtomicElement IRON;
  static const AtomicElement NICKEL;
  static const AtomicElement COBALT;
  static const AtomicElement COPPER;
  static const AtomicElement ZINC;
  static const AtomicElement GALLIUM;
  static const AtomicElement GERMANIUM;
  static const AtomicElement ARSENIC;
  static const AtomicElement SELENIUM;
  static const AtomicElement BROMINE;
  static const AtomicElement KRYPTON;
  static const AtomicElement RUBIDIUM;
  static const AtomicElement STRONTIUM;
  static const AtomicElement YTTRIUM;
  static const AtomicElement ZIRCONIUM;
  static const AtomicElement NIOBIUM;
  static const AtomicElement MOLYBDENUM;
  static const AtomicElement TECHNETIUM;
  static const AtomicElement RUTHENIUM;
  static const AtomicElement RHODIUM;
  static const AtomicElement PALLADIUM;
  static const AtomicElement SILVER;
  static const AtomicElement CADMIUM;
  static const AtomicElement INDIUM;
  static const AtomicElement TIN;
  static const AtomicElement ANTIMONY;
  static const AtomicElement TELLURIUM;
  static const AtomicElement IODINE;
  static const AtomicElement XENON;
  static const AtomicElement CAESIUM;
  static const AtomicElement BARIUM;
  static const AtomicElement LANTHANUM;
  static const AtomicElement CERIUM;
  static const AtomicElement PRASEODYMIUM;
  static const AtomicElement NEODYMIUM;
  static const AtomicElement PROMETHIUM;
  static const AtomicElement SAMARIUM;
  static const AtomicElement EUROPIUM;
  static const AtomicElement GADOLINIUM;
  static const AtomicElement TERBIUM;
  static const AtomicElement DYSPROSIUM;
  static const AtomicElement HOLMIUM;
  static const AtomicElement ERBIUM;
  static const AtomicElement THULIUM;
  static const AtomicElement YTTERBIUM;
  static const AtomicElement LUTETIUM;
  static const AtomicElement HAFNIUM;
  static const AtomicElement TANTALUM;
  static const AtomicElement TUNGSTEN;
  static const AtomicElement RHENIUM;
  static const AtomicElement OSMIUM;
  static const AtomicElement IRIDIUM;
  static const AtomicElement PLATINUM;
  static const AtomicElement GOLD;
  static const AtomicElement MERCURY;
  static const AtomicElement THALLIUM;
  static const AtomicElement LEAD;
//  static const AtomicElement URANIUM;
  static const AtomicElement BISMUTH;
  
  const index2 z;
  const real mass;
  const std::string symbol;
  const std::string name;
  
  static const std::unordered_map<std::string,AtomicElement> elements_by_symbol;

  /// Returns atomic element for a given mass (in [u] units)
  static const AtomicElement & by_mass(const core::real mass);

  /** @brief  Returns atomic element for a given element symbol string.
   * This method looks in elements_by_symbol map. If not found, it tries other
   * lower/upper case spelling variants
   */
  static const AtomicElement & by_symbol(const std::string &symbol);

  friend std::ostream& operator<< (std::ostream &out, const AtomicElement &e);

  static std::unordered_map<std::string,AtomicElement> create_map() {

    std::unordered_map<std::string,AtomicElement> tmp;
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::DUMMY.symbol,AtomicElement::DUMMY));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::HYDROGEN.symbol,AtomicElement::HYDROGEN));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::DEUTERIUM.symbol,AtomicElement::DEUTERIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::HELIUM.symbol,AtomicElement::HELIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::LITHIUM.symbol,AtomicElement::LITHIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::BERYLLIUM.symbol,AtomicElement::BERYLLIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::BORON.symbol,AtomicElement::BORON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CARBON.symbol,AtomicElement::CARBON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::NITROGEN.symbol,AtomicElement::NITROGEN));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::OXYGEN.symbol,AtomicElement::OXYGEN));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::FLUORINE.symbol,AtomicElement::FLUORINE));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::NEON.symbol,AtomicElement::NEON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SODIUM.symbol,AtomicElement::SODIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::MAGNESIUM.symbol,AtomicElement::MAGNESIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ALUMINIUM.symbol,AtomicElement::ALUMINIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SILICON.symbol,AtomicElement::SILICON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::PHOSPHORUS.symbol,AtomicElement::PHOSPHORUS));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SULFUR.symbol,AtomicElement::SULFUR));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CHLORINE.symbol,AtomicElement::CHLORINE));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::POTASSIUM.symbol,AtomicElement::POTASSIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ARGON.symbol,AtomicElement::ARGON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CALCIUM.symbol,AtomicElement::CALCIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SCANDIUM.symbol,AtomicElement::SCANDIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TITANIUM.symbol,AtomicElement::TITANIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::VANADIUM.symbol,AtomicElement::VANADIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CHROMIUM.symbol,AtomicElement::CHROMIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::MANGANESE.symbol,AtomicElement::MANGANESE));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::IRON.symbol,AtomicElement::IRON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::NICKEL.symbol,AtomicElement::NICKEL));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::COBALT.symbol,AtomicElement::COBALT));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::COPPER.symbol,AtomicElement::COPPER));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ZINC.symbol,AtomicElement::ZINC));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::GALLIUM.symbol,AtomicElement::GALLIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::GERMANIUM.symbol,AtomicElement::GERMANIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ARSENIC.symbol,AtomicElement::ARSENIC));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SELENIUM.symbol,AtomicElement::SELENIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::BROMINE.symbol,AtomicElement::BROMINE));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::KRYPTON.symbol,AtomicElement::KRYPTON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::RUBIDIUM.symbol,AtomicElement::RUBIDIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::STRONTIUM.symbol,AtomicElement::STRONTIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::YTTRIUM.symbol,AtomicElement::YTTRIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ZIRCONIUM.symbol,AtomicElement::ZIRCONIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::NIOBIUM.symbol,AtomicElement::NIOBIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::MOLYBDENUM.symbol,AtomicElement::MOLYBDENUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TECHNETIUM.symbol,AtomicElement::TECHNETIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::RUTHENIUM.symbol,AtomicElement::RUTHENIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::RHODIUM.symbol,AtomicElement::RHODIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::PALLADIUM.symbol,AtomicElement::PALLADIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SILVER.symbol,AtomicElement::SILVER));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CADMIUM.symbol,AtomicElement::CADMIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::INDIUM.symbol,AtomicElement::INDIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TIN.symbol,AtomicElement::TIN));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ANTIMONY.symbol,AtomicElement::ANTIMONY));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TELLURIUM.symbol,AtomicElement::TELLURIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::IODINE.symbol,AtomicElement::IODINE));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::XENON.symbol,AtomicElement::XENON));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CAESIUM.symbol,AtomicElement::CAESIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::BARIUM.symbol,AtomicElement::BARIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::LANTHANUM.symbol,AtomicElement::LANTHANUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::CERIUM.symbol,AtomicElement::CERIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::PRASEODYMIUM.symbol,AtomicElement::PRASEODYMIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::NEODYMIUM.symbol,AtomicElement::NEODYMIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::PROMETHIUM.symbol,AtomicElement::PROMETHIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::SAMARIUM.symbol,AtomicElement::SAMARIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::EUROPIUM.symbol,AtomicElement::EUROPIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::GADOLINIUM.symbol,AtomicElement::GADOLINIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TERBIUM.symbol,AtomicElement::TERBIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::DYSPROSIUM.symbol,AtomicElement::DYSPROSIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::HOLMIUM.symbol,AtomicElement::HOLMIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::ERBIUM.symbol,AtomicElement::ERBIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::THULIUM.symbol,AtomicElement::THULIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::YTTERBIUM.symbol,AtomicElement::YTTERBIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::LUTETIUM.symbol,AtomicElement::LUTETIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::HAFNIUM.symbol,AtomicElement::HAFNIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TANTALUM.symbol,AtomicElement::TANTALUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::TUNGSTEN.symbol,AtomicElement::TUNGSTEN));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::RHENIUM.symbol,AtomicElement::RHENIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::OSMIUM.symbol,AtomicElement::OSMIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::IRIDIUM.symbol,AtomicElement::IRIDIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::PLATINUM.symbol,AtomicElement::PLATINUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::GOLD.symbol,AtomicElement::GOLD));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::MERCURY.symbol,AtomicElement::MERCURY));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::THALLIUM.symbol,AtomicElement::THALLIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::LEAD.symbol,AtomicElement::LEAD));
//    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::URANIUM.symbol,AtomicElement::URANIUM));
    tmp.insert(std::pair<std::string,AtomicElement>(AtomicElement::BISMUTH.symbol,AtomicElement::BISMUTH));
    return tmp;
  }

private:
  AtomicElement(const core::index2 z,const real mass,const std::string symbol,const std::string name) : z(z), mass(mass), symbol(symbol), name(name) {}
  static utils::Logger logger;
};


std::ostream& operator<< (std::ostream &out, const AtomicElement &e);


}
}

#endif
