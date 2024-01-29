/** \file io_utils.hh
 * @brief Provides several I/O utility methods.
 */

#ifndef UTILS_io_utils_HH
#define UTILS_io_utils_HH

#include <zlib.h>
#include <sys/stat.h>

#include <memory>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>

namespace utils {

/** @brief Character used to separate sub-directories in a file path
 */
#ifdef OS_WINDOWS
static const char dir_separator = '\\';
#else
static const char dir_separator = '/';
#endif

/** @brief Joins two parts of a path into a single one using the separator appropriate to the OS
 *
 * @param p1 - the first part of the path, e.g. a directory
 * @param p2 - the second part of the path, e.g. a file name
 * @return a full path
 */
std::string join_paths(const std::string& p1, const std::string& p2);

/** @brief Lists the contents of a directory.
 *
 * This method actually wraps the glob() function from glibc
 * @param mask - a relative path with a file mask, e.g. <code>./some/dir/1*.pdb</code>
 * @return a vector holding file names
 * @see glibc site http://www.gnu.org/software/libc/manual/html_node/Globbing.html
 */
std::vector<std::string> glob(const std::string & mask);

/** @brief Returns a smart pointer to an output stream.
 *
 * If the given file name is an empty string, simply a pointer to std::cout is returned.
 * Otherwise the method returns a pointer to a newly opened file. Alternatively,
 * user may explicitly ask to write to std::cout or to std::cerr
 * by passing "stdout" or "stderr" as a file name.
 *
 * @param fname - name of the output file
 * @return an output stream
 */
std::shared_ptr<std::ostream> out_stream(const std::string &fname);

/** @brief Associates a disk file with a given input stream if the file exists.
 *
 * Otherwise, an exception is thrown. This method should be used as follows:
 * <code>
  std::ifstream in;
  in_stream(fname,in); // throws an exception if anything goes wrong
 *
 * </code>
 *
 * @param fname - name of the input file
 * @param in_stream - the stream  the file will be associated with
 * @param mode - flags describing the requested i/o mode for the file.
 */
void in_stream(const std::string &fname, std::ifstream & in_stream, std::ios_base::openmode mode = std::ios_base::in);

/** @brief Loads the whole content of a text file into a string.
 *
 * @param fname - name of the input file
 * @return a string containing all the characters from the given file
 */
std::string load_text_file(const std::string &fileName);

/** @brief Loads the whole content of a text file into the given string.
 *
 * @param fname - name of the input file
 * @param sink  - where the data should go; the previous content of this string is lost
 * @return a string containing all the characters from the given file
 */
void load_text_file(const std::string &fileName, std::string &sink);

/** @brief Loads properties file (the standard JAVA file format)
 *
 * @param fname - name of the input file
 * @param replace_underscores_with_spaces - if true, any '_' character will be replaced with ' '
 * @return a map that binds each string key to a vector of string tokens
 */
std::map<std::string,std::vector<std::string>> read_properties_file(const std::string& fname, const bool replace_underscores_with_spaces = false);

/** @brief Loads properties file (the standard JAVA file format).
 *
 * This method creates a new map and fills it with the data found in a given file.
 * The simple example loads a file and prints its content on a screen:
 * \include ex_read_properties_file.cc
 * @param fname - name of the input file
 * @param storage_map - container where the data will be stored
 * @param replace_undersores_with_spaces - if true, any '_' character will be replaced with ' '
 * @return a map that binds each string key to a vector of string tokens
 */
std::map<std::string, std::vector<std::string>> & read_properties_file(const std::string& fname,
    std::map<std::string, std::vector<std::string>> & storage_map, const bool replace_undersores_with_spaces);


std::string & load_binary_file(const std::string & filename, std::string & buffer);

std::string & unzip_string(const std::string& str,std::string & dest);
std::string & zip_string(const std::string& str, std::string& dest, int compressionlevel = Z_BEST_COMPRESSION);
std::stringstream & unzip_string(const std::string& str,std::stringstream & dest);

std::string & ungzip_string( const std::string& compressedBytes, std::string& uncompressedBytes );
std::stringstream & ungzip_string( const std::string& compressedBytes, std::stringstream& uncompressedBytes );

/** @brief Reads a file with strings and returns all of them as a vector instance.
 *
 * This method is intended to read a listfile i.e. a file providing a list of files
 * @param fname - an input file name
 * @return a vector holding file names
 */
std::vector<std::string> read_listfile(const std::string & fname);

/** @brief Reads a file with strings and places the strings in a given vector.
 *
 * This method is intended to read a listfile i.e. a file providing a list of files
 * @param fname - an input file name
 * @param destination - where to store strings read from a list-file
 */
void read_listfile(const std::string & fname,std::vector<std::string> & destination);

/** @brief Tests if a given file exists.
 *
 * @param fname - name of the input file
 * @return true if the requested file exists
 */
bool if_file_exists(const std::string& fname);

/// Returns a file name extracted from a path
std::string basename(const std::string& str);

/// Returns a file name root and extension (just the last extension is extracted from the root name)
inline std::pair<std::string, std::string> root_extension(const std::string& str) {

  return std::make_pair(str.substr(0, str.find_last_of(".")), str.substr(str.find_last_of(".") + 1));
}

/** @brief Removes extension from a file name
 *
 * @param fname - a file name to work on
 * @param extensions - a vector of extensions to be removed from the file name. Each of the extensions is tested
 * in the order as they appear in the vector.
 */
void trim_extensions(std::string & fname ,std::vector<std::string> extensions = {".pdb.gz",".pdb",".gz",".fasta",".pir"});

/// Returns a directory name extracted from a path
inline std::string pathname(const std::string& str) {

  return str.substr(0, str.find_last_of("/\\"));
}

/** @brief Tries to find a file in certain locations and returns true when succeeded.
 *
 * This method tests:
 *    - the given file name in a local directory
 *    - the given file name in the given directory
 *    - basename of the given file name in the given directory
 * The actual name of the file that has been found will be stored in the <code>result</code> variable
 * @param fname - name of the file to look for
 * @param result - if found, resulting file name (possibly with the correct path) will be stored at this reference
 * @param path - database path to check
 * @return true if the requested file exists
 */
bool find_file(const std::string & fname, std::string & result, const std::string & path);

/** @brief Creates a nicely looking timestamp.
 *
 * Returns a string with a timestamp, e.g. "26-11-2014,0:6:17"
 */
std::string time_stamp();

}
/**
 * \example ex_read_properties_file.cc
 */
#endif

