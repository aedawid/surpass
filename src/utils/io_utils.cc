#include <glob.h>
#include <zlib.h>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <sys/stat.h>

#include <memory>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <utils/io_utils.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

namespace utils {

utils::Logger io_utils_logger("io_utils");

bool if_file_exists(const std::string& fname) {
  struct stat buffer;
  return (stat(fname.c_str(), &buffer) == 0);
}

std::string basename(const std::string& str) {

  return str.substr(str.find_last_of("/\\") + 1);
}

std::string join_paths(const std::string& p1, const std::string& p2) {

  std::string tmp = p1;
  if (p1.back() != dir_separator) {
    tmp += dir_separator;
    return (tmp + p2);
  } else return (p1 + p2);
}

std::vector<std::string> glob(const std::string & mask) {

  io_utils_logger << utils::LogLevel::FILE << "listing files: " << mask << "\n";

  std::vector<std::string> out;
  glob_t globbuf;
  int err = ::glob(mask.c_str(), 0, NULL, &globbuf);
  if (err == 0) {
    std::cerr << globbuf.gl_pathc<<" files found\n";
    for (size_t i = 0; i < globbuf.gl_pathc; i++) {
      out.push_back(globbuf.gl_pathv[i]);
    }

    globfree(&globbuf);
  }
  return out;
}

std::vector<std::string> read_listfile(const std::string & fname) {

  std::vector<std::string> out;
  read_listfile(fname, out);

  return out;
}

void read_listfile(const std::string & fname, std::vector<std::string> & destination) {

  io_utils_logger.log(utils::LogLevel::FILE, "Reading a list-file named: ", fname, "\n");
  std::ifstream in;
  in_stream(fname, in);
  size_t cnt = 0;
  std::string line, token;
  try {
    while (std::getline(in, line)) {
      if (line.length() < 2) continue;
      if (line[0] == '#') continue;
      std::stringstream ss(line);
      ss >> token;
      destination.push_back(token);
      cnt++;
    }
  } catch (std::ifstream::failure & e) {
    if ((!in.fail()) || (!in.eof()) || (in.bad()))
      io_utils_logger.log(utils::LogLevel::CRITICAL, "Exception captured: ", e.what(), "\n");
  }
  io_utils_logger.log(utils::LogLevel::FILE, "found ", cnt, " file names in ", fname, "\n");
}

std::shared_ptr<std::ostream> out_stream(const std::string &fname) {

  if (fname.length() > 0) {
    if ((fname.compare("stderr") == 0) || (fname.compare("cerr") == 0)) {
      std::shared_ptr<std::ostream> p(&std::cerr, [](std::ostream*){});
      return p;
    }
    if ((fname.compare("stdout") == 0) || (fname.compare("cout") == 0)) {
      std::shared_ptr<std::ostream> p(&std::cout, [](std::ostream*){});
      return p;
    }
    if ((fname.compare("null") == 0) || (fname.compare("/dev/null") == 0)) {
      std::shared_ptr<std::ostream> p(new std::ofstream(0));
    }
    std::shared_ptr<std::ostream> p(new std::ofstream(fname));
    io_utils_logger.log(utils::LogLevel::FILE, "opening a file ", fname, "\n");
    return p;
  } else {
    std::shared_ptr<std::ostream> p(new std::ostream(std::cout.rdbuf()));
    return p;
  }
}

void in_stream(const std::string &fname, std::ifstream & in_stream, std::ios_base::openmode mode) {

  in_stream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  if (fname.length() > 0) if (if_file_exists(fname)) {
    in_stream.open(fname, mode);
    if (in_stream.is_open()) return;
  }
  io_utils_logger.log(utils::LogLevel::FILE, "Can't find file: ", fname, "\n");
}

std::string load_text_file(const std::string & fileName) {

  io_utils_logger.log(utils::LogLevel::FILE, "Attempting to open a text file: ", fileName, "\n");
  if (!if_file_exists(fileName)) {
    io_utils_logger.log(utils::LogLevel::SEVERE, "Can't find file: ", fileName, ", returning null-pointer!\n\n");
    throw std::runtime_error("Can't find file: " + fileName + ", returning null-pointer!\n");
  }

  std::ifstream ifs(fileName.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (!ifs) {
    io_utils_logger.log(utils::LogLevel::FILE, "Can't find file: " , fileName , ", returning null-pointer!\n\n");
    throw std::runtime_error("Can't find file: " + fileName + ", returning null-pointer!\n");
  }
  std::ifstream::pos_type fileSize = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  std::vector<char> bytes(fileSize);
  ifs.read(&bytes[0], fileSize);

  return std::string(&bytes[0], fileSize);
}

void load_text_file(const std::string &fileName, std::string &sink) {

  io_utils_logger.log(utils::LogLevel::FILE, "Attempting to open a text file: ", fileName, "\n");

  if (!if_file_exists(fileName)) {
    io_utils_logger.log(utils::LogLevel::SEVERE, "Can't find file: ", fileName, "returning null-pointer!\n\n");
    throw std::runtime_error("Can't find file: " + fileName + ", returning null-pointer!\n");
  }

  std::ifstream ifs(fileName.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (!ifs) {
    io_utils_logger.log(utils::LogLevel::SEVERE, "Can't find file: " , fileName , ", returning null-pointer!\n\n");
    throw std::runtime_error("Can't find file: " + fileName + ", returning null-pointer!\n");
  }
  std::ifstream::pos_type fileSize = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  std::vector<char> bytes(fileSize);
  ifs.read(&bytes[0], fileSize);
  sink.assign(bytes.begin(), bytes.end());
}

std::map<std::string, std::vector<std::string>> & read_properties_file(const std::string& fname,
    std::map<std::string, std::vector<std::string>> & storage_map, const bool replace_undersores_with_spaces) {

  std::ifstream in;
  in_stream(fname, in);
  size_t cnt = 0;
  std::string line, token;
  try {
    while (std::getline(in, line)) {
      if (line.length() < 2) continue;
      if (line[0] == '#') continue;
      if (line[0] == '!') continue;
      std::vector<std::string> tokens = split(line, ':');
      if (tokens.size() != 2) continue;
      std::vector<std::string> tokens2 = split(tokens[1], ' ');
      if (replace_undersores_with_spaces) std::replace(tokens[0].begin(), tokens[0].end(), '_', ' ');
      if (replace_undersores_with_spaces)
        for (std::string &s : tokens2)
          std::replace(s.begin(), s.end(), '_', ' ');
      storage_map[tokens[0]] = tokens2;
      cnt++;
    }
  } catch (std::ifstream::failure &e) {
    if ((!in.fail()) || (!in.eof()) || (in.bad()))
      io_utils_logger.log(utils::LogLevel::CRITICAL, "Exception captured: ", e.what(), "\n");
  }
  io_utils_logger.log(utils::LogLevel::FILE, "found ", cnt, " tokens in ", fname, "\n");

  return storage_map;
}

std::map<std::string, std::vector<std::string>> read_properties_file(const std::string& fname,
    const bool replace_underscores_with_spaces) {

  std::map<std::string, std::vector<std::string>> output;
  return read_properties_file(fname, output, replace_underscores_with_spaces);
}

std::string & load_binary_file(const std::string & filename, std::string & buffer) {

  if (!if_file_exists(filename)) {
    io_utils_logger.log(utils::LogLevel::SEVERE, "Can't find file: ", filename, ", returning null-pointer!\n\n");
    throw std::runtime_error("Can't find file: " + filename + ", returning null-pointer!\n");
  }
  std::ifstream input(filename, std::ios::binary);
  std::copy(std::istreambuf_iterator<char>(input), (std::istreambuf_iterator<char>()),
      std::inserter(buffer, buffer.begin()));
  return buffer;
}

std::string & zip_string(const std::string& str, std::string& dest, int compressionlevel) {

  z_stream zs;                        // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (deflateInit(&zs, compressionlevel) != Z_OK) throw(std::runtime_error("deflateInit failed while compressing."));

  zs.next_in = (Bytef*) str.data();
  zs.avail_in = str.size();           // set the z_stream's input

  int ret;
  char outbuffer[32768];
  dest.clear();

  // retrieve the compressed bytes blockwise
  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = deflate(&zs, Z_FINISH);

    if (dest.size() < zs.total_out) {
      // append the block to the output string
      dest.append(outbuffer, zs.total_out - dest.size());
    }
  } while (ret == Z_OK);

  deflateEnd(&zs);

  if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
    std::ostringstream oss;
    oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }

  return dest;
}

std::string & unzip_string(const std::string& str, std::string & dest) {

  z_stream zs;                        // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (inflateInit(&zs) != Z_OK) throw(std::runtime_error("inflateInit failed while decompressing."));

  zs.next_in = (Bytef*) str.data();
  zs.avail_in = str.size();

  int ret;
  char outbuffer[32768];
  dest.clear();

  // get the decompressed bytes blockwise using repeated calls to inflate
  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);
    ret = inflate(&zs, 0);
    if (dest.size() < zs.total_out) dest.append(outbuffer, zs.total_out - dest.size());

  } while (ret == Z_OK);

  inflateEnd(&zs);

  if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }

  return dest;
}

std::string& ungzip_string(const std::string& compressedBytes, std::string& uncompressedBytes) {

  if (compressedBytes.size() == 0) return uncompressedBytes;

  uncompressedBytes.clear();

  unsigned full_length = compressedBytes.size();
  unsigned half_length = compressedBytes.size() / 2;

  unsigned uncompLength = full_length;
  char* uncomp = (char*) calloc(sizeof(char), uncompLength);

  z_stream strm;
  strm.next_in = (Bytef *) compressedBytes.c_str();
  strm.avail_in = compressedBytes.size();
  strm.total_out = 0;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;

  bool done = false;
  int ret = 0;
  if ((ret = inflateInit2(&strm, (16+MAX_WBITS))) != Z_OK) {
    free(uncomp);
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << strm.msg;
    throw(std::runtime_error(oss.str()));
  }

  while (!done) {
    // If our output buffer is too small
    if (strm.total_out >= uncompLength) {
      // Increase size of output buffer
      char* uncomp2 = (char*) calloc(sizeof(char), uncompLength + half_length);
      memcpy(uncomp2, uncomp, uncompLength);
      uncompLength += half_length;
      free(uncomp);
      uncomp = uncomp2;
    }

    strm.next_out = (Bytef *) (uncomp + strm.total_out);
    strm.avail_out = uncompLength - strm.total_out;

    // Inflate another chunk.
    int err = inflate(&strm, Z_SYNC_FLUSH);
    if (err == Z_STREAM_END) done = true;
    else if (err != Z_OK) {
      std::ostringstream oss;
      oss << "Exception during zlib decompression: (" << err << ") " << strm.msg;
      throw(std::runtime_error(oss.str()));
    }
  }

  if ((ret = inflateEnd(&strm)) != Z_OK) {
    free(uncomp);
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << strm.msg;
    throw(std::runtime_error(oss.str()));
  }

  for (size_t i = 0; i < strm.total_out; ++i) {
    uncompressedBytes += uncomp[i];
  }
  free(uncomp);
  return uncompressedBytes;
}

std::stringstream & unzip_string(const std::string& str, std::stringstream & dest) {

  z_stream zs;                        // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (inflateInit(&zs) != Z_OK) throw(std::runtime_error("inflateInit failed while decompressing."));

  zs.next_in = (Bytef*) str.data();
  zs.avail_in = str.size();

  int ret;
  char outbuffer[32768];
  dest.clear();

  // get the decompressed bytes block-wise using repeated calls to inflate
  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = inflate(&zs, 0);
    dest << outbuffer;

  } while (ret == Z_OK);

  inflateEnd(&zs);

  if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }

  return dest;
}

std::stringstream & ungzip_string(const std::string& compressedBytes, std::stringstream & uncompressedBytes) {

  if (compressedBytes.size() == 0) return uncompressedBytes;

  uncompressedBytes.clear();

  unsigned full_length = compressedBytes.size();
  unsigned half_length = compressedBytes.size() / 2;

  unsigned uncompLength = full_length;
  char* uncomp = (char*) calloc(sizeof(char), uncompLength);

  z_stream strm;
  strm.next_in = (Bytef *) compressedBytes.c_str();
  strm.avail_in = compressedBytes.size();
  strm.total_out = 0;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;

  bool done = false;
  int ret = 0;
  if ((ret = inflateInit2(&strm, (16+MAX_WBITS))) != Z_OK) {
    free(uncomp);
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << strm.msg;
    throw(std::runtime_error(oss.str()));
  }

  while (!done) {
    // If our output buffer is too small
    if (strm.total_out >= uncompLength) {
      // Increase size of output buffer
      char* uncomp2 = (char*) calloc(sizeof(char), uncompLength + half_length);
      memcpy(uncomp2, uncomp, uncompLength);
      uncompLength += half_length;
      free(uncomp);
      uncomp = uncomp2;
    }

    strm.next_out = (Bytef *) (uncomp + strm.total_out);
    strm.avail_out = uncompLength - strm.total_out;

    // Inflate another chunk.
    int err = inflate(&strm, Z_SYNC_FLUSH);
    if (err == Z_STREAM_END) done = true;
    else if (err != Z_OK) {
      std::ostringstream oss;
      oss << "Exception during zlib decompression: (" << err << ") " << strm.msg;
      throw(std::runtime_error(oss.str()));
    }
  }

  if ((ret = inflateEnd(&strm)) != Z_OK) {
    free(uncomp);
    std::ostringstream oss;
    oss << "Exception during zlib decompression: (" << ret << ") " << strm.msg;
    throw(std::runtime_error(oss.str()));
  }

//  for (size_t i = 0; i < strm.total_out; ++i) {
//    uncompressedBytes << uncomp[i];
//  }
  uncompressedBytes << uncomp;

  free(uncomp);
  return uncompressedBytes;
}

bool find_file(const std::string & fname, std::string & result, const std::string & path) {

  if (if_file_exists(fname)) {
    result = fname;
    return true;
  }
  std::string lastname = basename(fname);

  result = join_paths(path, fname);
  if (if_file_exists(result)) return true;

  result = join_paths(path, basename(fname));
  if (if_file_exists(result)) return true;

  return false;
}

std::string time_stamp() {

  time_t ltime;
  struct tm *Tm;
  ltime = time(NULL);
  Tm = localtime(&ltime);
  return utils::string_format("%d-%d-%d,%d:%d:%d", Tm->tm_mday, Tm->tm_mon + 1, Tm->tm_year + 1900, Tm->tm_hour,
      Tm->tm_min, Tm->tm_sec);
}

void trim_extensions(std::string & fname,std::vector<std::string> extensions ) {

  for (const std::string &ext : extensions)
    if (fname.find(ext) != std::string::npos) utils::replace_substring(fname, ext, "");
}

}

