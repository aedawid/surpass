/** @file Provides SequenceFilter type and several handy filter objects.
 */
#ifndef CORE_DATA_SEQUENCE_SequenceFilter_H
#define CORE_DATA_SEQUENCE_SequenceFilter_H

#include <string>
#include <functional>

#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/sequence_utils.hh>
#include <core/alignment/on_alignment_computations.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace sequence {

/** @brief SequenceFilter returns true if a given sequence satisfies certain criteria
 *
 */
class SequenceFilter : public std::unary_function<std::shared_ptr<Sequence>,bool> {
public:
  /// Mandatory virtual destructor
  virtual ~SequenceFilter() { }
  /// Apply a filter to the reference to a sequence
  virtual bool operator()(const Sequence & s) const = 0;
  /// Apply a filter to a pointer to a sequence
  virtual bool operator()(const std::shared_ptr<Sequence> s) const = 0;
  /// Apply to a sequence given as a string
  virtual bool operator()(const std::string & header, const std::string & sequence) const = 0;
};

/* @brief Returns true if a sequence header contains a certain substring
 *
 */
class FindInSequenceName : public SequenceFilter {
public:
  /// The filter will pass only these sequences whose name contains a given substring
  FindInSequenceName(const std::string & wanted_substr) { wanted_substrings.push_back(wanted_substr); }

  /* @brief The filter will pass only these sequences whose name contains at least one of the given substrings
   * @param  wanted_substr - a vector of possible substring that make a sequence eligible
   */
  FindInSequenceName(const std::vector<std::string> & wanted_substr) {
    std::copy(wanted_substr.begin(), wanted_substr.end(), wanted_substrings.begin());
  }

  virtual inline bool operator()(const Sequence & s) const {

    for (const std::string si : wanted_substrings)
      if (s.header().find(si)!=std::string::npos) return true;
    return false;
  }

  virtual inline bool operator()(const std::shared_ptr<Sequence> s) const {

    for (const std::string si : wanted_substrings)
      if (s->header().find(si)!=std::string::npos) return true;

    return false;
  }

  virtual inline bool operator()(const std::string & header, const std::string & sequence) const {

    for (const std::string si : wanted_substrings)
      if (header.find(si)!=std::string::npos) return true;
    return false;
  }

  /** @brief Loads the list of wanted substring (subsequences) from a file.
   *
   * @param file_name - the name of an input file
   */
  void substrings_from_file(const std::string & file_name);

  /** @brief Loads the list of wanted substring (subsequences) from a stream.
   *
   * @param file_name - the input stream
   */
  void substrings_from_stream(std::istream &in_stream);

private:
  std::vector<std::string> wanted_substrings;
  static utils::Logger logger;
};

/** @brief Returns true if a given sequence is a protein
 *
 */
class IsProteinSequence : public SequenceFilter {
public:

/** @brief Creates a filter that passes  protein sequences.
 * A protein sequence string is recognized by having at least <code>aa_fraction</code> fraction of amino acids.
 * By default this parameter is set to 0.9 so sequences taken from PDB chains pass the filter even if they contain ligand residues
 */
  IsProteinSequence(const real aa_fraction = 0.9) : aa_fraction_(aa_fraction) {}

  virtual inline bool operator()(const Sequence & s) const { return s.is_protein_sequence; }

  virtual inline bool operator()(const std::shared_ptr<Sequence> s) const { return operator()(*s); }

  virtual inline bool operator()(const std::string & header, const std::string & sequence) const {

    return count_aa(sequence) / float(sequence.length())> aa_fraction_;
  }

private:
  core::real aa_fraction_;
};

/** @brief Returns true if a given sequence is a nucleic acid sequence
 *
 */
class IsNucleicSequence : public SequenceFilter {
public:

  IsNucleicSequence() {}

  virtual inline bool operator()(const Sequence & s) const { return s.is_nucleic_sequence; }

  virtual inline bool operator()(const std::shared_ptr<Sequence> s) const {  return operator()(*s); }

  virtual inline bool operator()(const std::string & header, const std::string & sequence) const {

    return count_nt(sequence) == sequence.size();
  }
};

/** @brief Returns true if a given sequence is not shorter than a given length (amino acid or nucleic acid residues)
 */
class NoShorterThan : public SequenceFilter {
public:

/** @brief Creates a filter that filters out too short sequences
 */
  NoShorterThan(const index2 at_least_that_long = 30) : min_length_(at_least_that_long) {}

  virtual inline bool operator()(const Sequence & s) const { return s.length()>= min_length_; }

  virtual inline bool operator()(const std::shared_ptr<Sequence> s) const { return s->length()>= min_length_; }

  virtual inline bool operator()(const std::string & header, const std::string & sequence) const {

    return sequence.length()>= min_length_;
  }

private:
  core::index2 min_length_;
};

/** @brief Returns true if a given sequence contains at least a given number of residue types
 *
 */
class IsNotTooSimple : public SequenceFilter {
public:

/** @brief Creates a filter that filters out too short sequences
 */
  IsNotTooSimple(const index2 n_distinct_res = 5) : n_distinct_res_(n_distinct_res) {}

  virtual inline bool operator()(const Sequence & s) const {

    return core::data::sequence::count_aa_types(s)>= n_distinct_res_;
  }

  virtual inline bool operator()(const std::shared_ptr<Sequence> s) const {
    return operator()(*s);
  }

  virtual inline bool operator()(const std::string & header, const std::string & sequence) const {

    return core::data::sequence::count_aa_types(sequence)>= n_distinct_res_;
  }

private:
  core::index2 n_distinct_res_;
};

/** @brief Returns true if a given sequence contains not too many unknown residues :
 * 'X' in one-letter code or "UNK", "UNG", "UNL" in three-letters code
 *
 */
class ByUnknownRatio : public SequenceFilter {
public:

/** @brief Creates a filter that filters out too short sequences
 */
  ByUnknownRatio(const real known_fraction = 0.9) : known_fraction_(known_fraction) {}

  virtual bool operator()(const Sequence & s) const;

  virtual inline bool operator()(const std::shared_ptr<Sequence> s) const { return operator()(*s); }

  virtual bool operator()(const std::string & header, const std::string & sequence) const;

private:
  core::real known_fraction_;
};

}
}
}

#endif
