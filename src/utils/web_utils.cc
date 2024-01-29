#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <iostream>
#include <iomanip>

#include <utils/web_utils.hh>
#include <utils/WebServer.hh>

#include <utils/Logger.hh>
#include <utils/options/Option.hh>
#include <utils/options/OptionParser.hh>
#include <utils/string_utils.hh>


#define PORT 80
#define RCVBUFSIZE 4096

namespace utils {

Logger logger("web_utils");

unsigned long resolve_name(const std::string & name) {

  struct hostent *host;

  if ((host = gethostbyname(name.c_str())) == NULL) {
    logger << LogLevel::SEVERE << "Can't resolve address: " << name << " - gethostbyname() failed\n";
    return 0;
  }

  return *((unsigned long *) host->h_addr);
}

void split_url(const std::string & url, std::string & first, std::string & second) {

  size_t i = 0;
  if (url[5] == '/' && url[6] == '/') i = 7;

  first.clear();
  for (; url[i] != '/'; i++)
    first += url[i];

  second.clear();
  for (; i < url.length(); i++)
    second += url[i];
}

std::string encode_url(const std::string &value) {

  std::ostringstream escaped;
  escaped.fill('0');
  escaped << std::hex;

  for (std::string::const_iterator i = value.begin(), n = value.end(); i != n; ++i) {
    std::string::value_type c = (*i);
    if (isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~') {
      escaped << c;
    } else if (c == ' ') {
      escaped << '+';
    } else {
      escaped << '%' << std::setw(2) << ((int) c) << std::setw(0);
    }
  }

  return escaped.str();
}

std::string decode_url(const std::string & url) {

  std::ostringstream decoded;
  decoded.fill('0');
  char *end = NULL;
  std::string substr ("  ");

  for (std::string::const_iterator i = url.begin(), n = url.end(); i != n; ++i) {
    std::string::value_type c = (*i);
    if(c!='%') {
      if(c=='+') {
        decoded << ' ';
        continue;
      }
      decoded << c;
    }
    else {
      substr[0] = (*(++i));
      substr[1] = (*(++i));
      char ascii = strtoul(substr.c_str(), &end, 16);
      decoded << ascii;
    }
  }

  return decoded.str();
}

std::string load_url(const std::string & url) {

  int sock, bytes;
  char buffer[RCVBUFSIZE];
  struct sockaddr_in server;
  std::stringstream ret;
  std::string address;
  std::string full_fname;
  split_url(url, address, full_fname);
  size_t pos = full_fname.find_last_of("/");
  std::string fname = full_fname.substr(pos);

  logger<<LogLevel::INFO<<"downloading "<<full_fname<<" from "<<address<<"\n";

  std::string command = "GET " + full_fname + " HTTP/1.1\r\nHost: " + address + "\r\n\r\n";

  /* create reliable stream socket */
  if ((sock = socket( PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    logger << LogLevel::SEVERE << "socket() failed\n";
    return ret.str();
  }

  /* Construct the server address structure */
  memset(&server, 0, sizeof(server));
  server.sin_family = AF_INET;
  server.sin_addr.s_addr = resolve_name(address);
  server.sin_port = htons(PORT);

  /* establish a connection to the server */
  if (connect(sock, (struct sockaddr *) &server, sizeof(server)) < 0) {
    logger << LogLevel::SEVERE << "connect() failed\n";
    return ret.str();
  }

  if (send(sock, command.c_str(), command.size(), 0) != int(command.size())) {
    logger << LogLevel::SEVERE << "send() sent a different number of bytes than expected\n";
    return ret.str();
  }

  int init = 0;
  while ((bytes = recv(sock, buffer, sizeof(buffer), 0)) > 0) {

    char *p;
    if (init) {
      p = buffer;
    } else {
      if ((p = (char *) memmem(buffer, bytes, "\r\n\r\n", 4))) {
        p += 4;
        init = 1;
      } else if ((p = (char *) memmem(buffer, bytes, "\n\n", 2))) {
        p += 2;
        init = 1;
      } else {
        continue;
      }
    }
    bytes -= (p - buffer);
    ret << buffer;
  }
  close(sock);

  return ret.str();
}

void log_request(Request* req) {

  logger << utils::LogLevel::INFO << "accessed: " << req->method << " " << req->path << "\n";

  if (logger.is_logable(utils::LogLevel::FINE)) {
    logger << utils::LogLevel::FINE << "Headers:" << "\n";
    std::map<std::string, std::string>::iterator iter;
    for (iter = req->headers.begin(); iter != req->headers.end(); ++iter) {
      logger << utils::LogLevel::FINE << iter->first << " = " << iter->second << "\n";
    }

    std::cout << "Query:" << "\n";
    for (iter = req->query.begin(); iter != req->query.end(); ++iter) {
      logger << utils::LogLevel::FINE << iter->first << " = " << iter->second << "\n";
    }

    std::cout << "Cookies: " << req->cookies.size() << std::endl;
    for (iter = req->cookies.begin(); iter != req->cookies.end(); ++iter) {
      logger << utils::LogLevel::FINE << iter->first << " = " << iter->second << "\n";
    }
  }
}

bool get_query_value(Request* req, const std::string & query_key, std::string & value) {

  auto a = req->query.find(query_key);
  if (a != req->query.end()) {
    value = utils::decode_url(req->query.find(query_key)->second);
    logger << utils::LogLevel::FINE << "value >" << value << "< for the key >" << query_key << "\n";
    return true;
  }
  return false;
}

std::string get_query_value(Request* req, const std::string & query_key) {

  auto a = req->query.find(query_key);
  if (a != req->query.end()) {
    std::string value = utils::decode_url(req->query.find(query_key)->second);
    logger << utils::LogLevel::FINE << "value >" << value << "< for the key >" << query_key << "\n";
    return value;
  }
  return "";
}

bool inject_http_option(Request* req, const std::string & query_key, utils::options::Option & option) {

  auto a = req->query.find(query_key);
  if (a != req->query.end()) {
    std::string value = utils::decode_url(req->query.find(query_key)->second);
    logger << utils::LogLevel::FINE << "value >" << value << "< for the key >" << query_key << "< assigned to option "
        << option.name << "\n";
    utils::options::OptionParser::get().inject_option(option, value);
    return true;
  }
  return false;
}

void show_error(Response* res, const std::string & msg) { res->send(msg); }

}
