#ifndef SIMULATIONS_OBSERVERS_PymolObserver_HH
#define SIMULATIONS_OBSERVERS_PymolObserver_HH

#include<iostream>
#include<arpa/inet.h>
#include<unistd.h>
#include<sys/socket.h>
#include<sys/types.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include <iostream>
#include <fstream>

#include <core/data/basic/Vec3.hh>
#include <simulations/systems/CartesianAtomsSimple.hh>
#include <simulations/observers/cartesian/AbstractPdbObserver.hh>

namespace simulations {
namespace observers {
namespace cartesian {

/** @brief Sends a conformation to PyMol program.
 *
 * Inside PyMol you must start a UDP server that listens to a port where the data is sent.
 * @param observed_object - system whose coordinates will be stored in the file
 * @param pdb_format_source - biomolecular structure that corresponds to the system.
 * @param address - address of the UDP server that listens to the PDB data (usually this is "127.0.0.1")
 * @param port - port to listen to , by default 65000
 * @see AbstractPdbObserver<C>
 */
template<typename C>
class PymolObserver : public AbstractPdbObserver<C> {
public:

  PymolObserver(const systems::CartesianAtomsSimple <C> &observed_object,
                const core::data::structural::Structure &pdb_format_source, const std::string address =
  "127.0.0.1", const size_t port = 65000) :
    AbstractPdbObserver<C>(observed_object, pdb_format_source) {

    sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    serv.sin_family = AF_INET;
    serv.sin_port = htons(65000);
    serv.sin_addr.s_addr = inet_addr(address.c_str());
    m = sizeof(serv);
  }

  virtual ~PymolObserver() {}

  /// Sends conformation to pymol as PDB formatted text
  virtual bool observe();

  /// Send some text to the server (server must be able to understand it!)
  bool observe(const std::string text);

  void bond(const core::index4 i_atom, const core::index4 j_atom) {
    send_line(utils::string_format("BOND bond id#%d, id#%d", i_atom, j_atom));
    sendto(sockfd, "Q", 1, 0, (struct sockaddr *) &serv, m);
  }

  /// The virtual method does nothing in this class
  virtual void finalize() {}

  virtual std::shared_ptr<std::ostream> output_stream() { return nullptr; }

  virtual void output_stream(std::shared_ptr<std::ostream> out) {}

private:
  int sockfd;
  struct sockaddr_in serv;
  socklen_t m = 0;
  core::index4  cnt = 0;

  void send_line(const std::string &line) {

    sendto(sockfd, line.c_str(), line.length(), 0, (struct sockaddr *) &serv, m);
    if (line.back() != '\n') sendto(sockfd, "\n", 1, 0, (struct sockaddr *) &serv, m);
  }
};

template<typename C>
bool PymolObserver<C>::observe() {

  ++cnt;

  std::stringstream oss;
  AbstractPdbObserver<C>::observed_object.write_pdb(oss, AbstractPdbObserver<C>::format_lines, cnt);
  while (!oss.eof()) {
    std::string line;
    getline(oss, line);
    send_line(line);
  }
  sendto(sockfd, "Q", 1, 0, (struct sockaddr *) &serv, m);
  return true;
}

template<typename C>
bool PymolObserver<C>::observe(const std::string text) {

  std::stringstream iss(text);
  while (!iss.eof()) {
    std::string line;
    getline(iss, line);
    send_line(line);
  }
  sendto(sockfd, "Q", 1, 0, (struct sockaddr *) &serv, m);
  return true;
}

}
}
}

#endif
