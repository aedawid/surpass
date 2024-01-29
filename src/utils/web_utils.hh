/** \file web_utils.hh
 * @brief Provides utilities for  web apps
 */

#ifndef UTILS_web_utils_HH
#define UTILS_web_utils_HH

#include <string>

#include <utils/WebServer.hh>

#include <utils/Logger.hh>
#include <utils/string_utils.hh>
#include <utils/options/Option.hh>

namespace utils {

unsigned long resolve_name(const std::string & name);

void split_url(const std::string & url, std::string & first, std::string & second);

std::string load_url(const std::string & url);

/** \brief Decodes a URL string that has been sent in an encoded form
 *
 * @param url - encoded URL
 * @return decoded URL string
 */
std::string decode_url(const std::string & url);

/** \brief Encodes a URL string by replacing non-ASCI characters with their codes.
 *
 * @param url - a URL to be encoded
 * @return encoded URL string
 */
std::string encode_url(const std::string & url);

/// Writes a log message about a given http request to the logger
void log_request(Request* req);

/** \brief Reads a value associated with a given key in a HTTP request.
 *
 * This method attempts to find a value for a key in a GET or POST request. For example, if a
 * GET request looked like: <code>http://www.site.com?protein_id=2gb1</code>, then
 * <code>
 *   get_query_value(req,"protein_id");
 * </code>
 * should return <code>2gb1</code> string.
 *
 * @param req - incoming http request
 * @param query_key - the key
 * @param value - placeholder to store the value
 * @return true if the given key was found in the request; false otherwise
 */
bool get_query_value(Request* req, const std::string & query_key, std::string & value);

/** \brief Reads a value associated with a given key in a HTTP request.
 *
 * This method attempts to find a value for a key in a GET or POST request. For example, if a
 * GET request looked like: <code>http://www.site.com?protein_id=2gb1</code>, then
 * <code>
 *   get_query_value(req,"protein_id");
 * </code>
 * should return <code>2gb1</code> string.
 *
 * @param req - incoming http request
 * @param query_key - the key
 * @return value associated with a given key
 */
std::string get_query_value(Request* req, const std::string & query_key);

/** \brief Registers an option=value pair in the OptionParser instance as it would be used at command-line.
 *
 * This method extracts from the given GET/POST query a value associated with the given string key. The value
 * is subsequently assigned to a command-line option.
 * @param req - incoming http request
 * @param query_key - GET/POST parameter identifier.
 * @param option - option to which a given value will be assigned
 * @return true if the given <code>query_key</code> has been found in the  given GET/POST query; false otherwise
 * (which means the option was not placed at the command line)
 */
bool inject_http_option(Request* req, const std::string & query_key, utils::options::Option & option);

/** \brief  Loads the given error message to a http package body.
 *
 * @param res - package that will be returned
 * @param msg - the error message text
 */
void show_error(Response* res, const std::string & msg);

}

#endif
