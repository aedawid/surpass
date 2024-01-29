#include <cmath>
#include <string>

#include <core/real.hh>
#include <core/index.hh>
#include <utils/string_utils.hh>

namespace utils {

#define BUFF_SIZE 4096
char sprintf_buffer[BUFF_SIZE];

const std::string letters("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

std::string string_format(const std::string & fmt_str, ...) {

  va_list ap;
  va_start(ap, fmt_str);
  vsnprintf(sprintf_buffer, BUFF_SIZE, fmt_str.c_str(), ap);
  va_end(ap);

  return std::string(sprintf_buffer);
}

std::string string_format(const char* fmt_str, ...) {

  va_list ap;
  va_start(ap, fmt_str);
  vsnprintf(sprintf_buffer, BUFF_SIZE, fmt_str, ap);
  va_end(ap);

  return std::string(sprintf_buffer);
}

std::string string_format(const char* fmt_str, std::vector<core::real> data) {

  std::stringstream ss;
  for (core::index2 i = 0; i < data.size(); ++i)
    ss << string_format(fmt_str, data[i]);

  return ss.str();
}

std::string& trim(std::string& str, const std::string delim) {

  std::string::size_type pos = str.find_last_not_of(delim);
  if (pos != std::string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(delim);
    if (pos != std::string::npos) str.erase(0, pos);
  } else str.erase(str.begin(), str.end());

  return str;
}

std::vector<std::string> split(const std::string &s, const char delim) {

  std::vector<std::string> tokens;
  split(s, tokens, delim);

  return tokens;
}

template<>
std::vector<std::string> & split<std::string>(const std::string &s, std::vector<std::string> &tokens,
    const char delim) {

  std::string s_copy(s);
  trim(s_copy);
  std::stringstream ss(s_copy);
  std::string item;
  while (std::getline(ss, item, delim)) {
    if (trim(item).size() == 0) continue;
    tokens.push_back(item);
  }

  return tokens;
}

std::vector<std::string> & split(const std::string &s, std::vector<std::string> &tokens, const core::index2 chunk_length) {

  for (unsigned i = 0; i < s.length(); i += chunk_length)
    tokens.push_back(s.substr(i, chunk_length));

  return tokens;
}

std::string& replace_substring(std::string& subject, const std::string& search, const std::string& replace) {

  if ((subject.length() == 0) || (search.length() == 0)) return subject;

  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != std::string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }

  return subject;
}

bool ends_with(const std::string & a, const std::string& b) {
  if (b.size() > a.size()) return false;
  return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}

void to_lower(std::string &s) {

  std::string::iterator i = s.begin();
  std::string::iterator end = s.end();

  while (i != end) {
    *i = std::tolower((unsigned char) *i);
    ++i;
  }
}

void to_upper(std::string &s) {

  std::string::iterator i = s.begin();
  std::string::iterator end = s.end();

  while (i != end) {
    *i = std::toupper((unsigned char) *i);
    ++i;
  }
}

double to_double(const char *p) {
  double r = 0.0;
  bool neg = false;
  while (*p == ' ')
    ++p;
  if (*p == '-') {
    neg = true;
    ++p;
  }
  while (*p >= '0' && *p <= '9') {
    r = (r * 10.0) + (*p - '0');
    ++p;
  }
  if (*p == '.') {
    double f = 0.0;
    int n = 0;
    ++p;
    while (*p >= '0' && *p <= '9') {
      f = (f * 10.0) + (*p - '0');
      ++p;
      ++n;
    }
    r += f / std::pow(10.0, n);
  }
  if (neg) {
    r = -r;
  }
  return r;
}

int to_int(const char *p) {

  int r = 0;
  bool neg = false;
  while (*p == ' ')
    ++p;
  if (*p == '-') {
    neg = true;
    ++p;
  }
  while (*p >= '0' && *p <= '9') {
    r = core::real((r * 10.0) + (*p - '0'));
    ++p;
  }

  if (neg) {
    r = -r;
  }
  return r;
}

bool is_integer(const std::string & s) {

  if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;
  char * p ;
  strtol(s.c_str(), &p, 10) ;

  return (*p == 0) ;
}

std::string format_paragraph(const std::vector<std::string> & words,
    const std::string & paragraph_pad, const std::string & line_pad, const core::index2 max_line_width) {

  std::string out = paragraph_pad;
  core::index2 lineLength = paragraph_pad.length();
  size_t pos = 0;
  for (const std::string & wrd : words) {
    std::string w(wrd);
    w = utils::trim(w);
    if ((pos = w.find("%N")) != w.npos) {
      w[pos] = ' ';
      w[pos + 1] = '\n';
      lineLength = 0;
    }
    if ((pos = w.find("%T")) != w.npos) {
      w[pos] = ' ';
      w[pos + 1] = '\t';
      lineLength += 4;
    }
    if (lineLength + w.length() + 1 > max_line_width) {
      // won't fit. Start a new line.
      if (lineLength != 0) {
        out += '\n';
        out += line_pad;
        lineLength = line_pad.length();
      }
      // no lead space
    } else {
      /* will fit */
      if (lineLength != 0) {
        // add lead space
        out += ' ';
        lineLength++;
      }
    }
    out += w;
    lineLength += w.length();
  } // end for

  return out;
}

bool has_suffix(const std::string &str, const std::string &suffix) {
  return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string format_paragraph(const std::vector<std::string> & words,
    const std::string & paragraph_pad, const std::string & line_pad, const int max_line_width, const int max_first_line_width) {

  std::string out = paragraph_pad;
  int lineLength = paragraph_pad.length();
  size_t pos = 0;
  bool is_first_line = true;
  for (const std::string & wrd : words) {
    std::string w(wrd);
    w = utils::trim(w);
    if ((pos = w.find("%N")) != w.npos) {
      w[pos] = ' ';
      w[pos + 1] = '\n';
      lineLength = 0;
    }
    if ((pos = w.find("%T")) != w.npos) {
      w[pos] = ' ';
      w[pos + 1] = '\t';
      lineLength += 4;
    }
    const size_t mx = (is_first_line) ? max_first_line_width : max_line_width;
    if (lineLength + w.length() + 1 > mx) {
      // won't fit. Start a new line.
      if (lineLength != 0) {
        out += '\n';
        out += line_pad;
        lineLength = line_pad.length();
        is_first_line = false;
      }
      // no lead space
    } else {
      /* will fit */
      if (lineLength != 0) {
        // add lead space
        out += ' ';
        lineLength++;
      }
    }
    out += w;
    lineLength += w.length();
  } // end for

  return out;
}

void random_alphanumeric(std::string & result) {

  static const char alphanum[] = "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";

  for (size_t i = 0; i < result.size(); ++i) result[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
}

std::pair<size_t,size_t> matching_bracket(const std::string & text, const size_t starting_pos, const std::string & bracket_type, const bool eat_brackets) {

  size_t beg_pos = text.find(bracket_type[0],starting_pos);
  short int state = 1; // one opened
  for (size_t i = beg_pos + 1; i < text.size(); ++i) {
    if (text[i] == bracket_type[0]) ++state;
    else if (text[i] == bracket_type[1]) --state;
    if(state == 0) { // stopping criteria : found a matching pair
      if (eat_brackets) return std::make_pair(beg_pos + 1, i);
      else return std::make_pair(beg_pos, i+1);
    }
  }

  return std::make_pair(0,0);
}

std::string to_base64(const std::string &in) {

  std::string out;

  int val=0, valb=-6;
  for (unsigned char c : in) {
    val = (val<<8) + c;
    valb += 8;
    while (valb>=0) {
      out.push_back("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[(val>>valb)&0x3F]);
      valb-=6;
    }
  }
  if (valb>-6) out.push_back("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[((val<<8)>>(valb+8))&0x3F]);
  while (out.size()%4) out.push_back('=');
  return out;
}

std::string from_base64(const std::string &in) {

  std::string out;

  std::vector<int> T(256,-1);
  for (int i=0; i<64; i++) T["ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"[i]] = i;

  int val=0, valb=-8;
  for (unsigned char c : in) {
    if (T[c] == -1) break;
    val = (val<<6) + T[c];
    valb += 6;
    if (valb>=0) {
      out.push_back(char((val>>valb)&0xFF));
      valb-=8;
    }
  }
  return out;
}

} // ~utils

