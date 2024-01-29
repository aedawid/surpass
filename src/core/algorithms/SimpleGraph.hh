#ifndef CORE_ALGORITHMS_SimpleGraph_HH
#define CORE_ALGORITHMS_SimpleGraph_HH

#include <vector>
#include <map>

namespace core {
namespace algorithms {

/** \brief Represents a simple graph.
 * <p>
 * A graph is a collection of points and lines connecting some (possibly empty) subset of them.
 * The points of a graph are referred as graph vertices.
 * Similarly, the lines connecting the vertices of a graph are most called edges.
 * A simple graph is an unweighted, undirected graph containing no graph loops or multiple edges.
 * </p>
 * <p>
 * This class represents a graph for which:</p>
 * <ul>
 * <li>  multiple edges are not allowed between vertices \f$\Longleftrightarrow \f$ (i.e. it is not a multigraph) </li>
 * <li>  edges do not support arrows \f$\Longleftrightarrow \f$ undirected graph </li>
 * <li>  Edge objects (of a generic type E) and vertices (of a generic type V) may be considered as
 * edge and vertex labels, respectively \f$\Longleftrightarrow \f$ labeled graph </li>
 * </ul>
 * </p>
 * @tparam <V> - A type for vertices
 * @tparam <E> - A type for edges
 * @see http://mathworld.wolfram.com/Graph.html for graph terminology and a brief description
 */
template<class V, class E>
class SimpleGraph {
public:

  /// a sole constructor creates an empty graph
  SimpleGraph() : logger("SimpleGraph") { }

  /** \brief Returns true if a given vertex belongs to this graph.
   * @param possible_vertex - a object that might be a vertex of this graph
   * @return true if possibleVertex really belongs to this graph
   */
  bool has_vertex(const V & possible_vertex) const { return (edgesForVertex.find(possible_vertex) != edgesForVertex.end()); }

  /**  \brief Adds a vertex into this graph.
   *
   * The method does not check if the graph already has this vertex! If you are not sure,
   * use hasVertex() to check if you really  need to insert.
   * @param vertex a new vertex to be inserted
   * @return integer index referring to the inserted vertex; the index may be further used to declare bonds within this molecule
   */
  core::index4 add_vertex(V & vertex) {

    vertices.push_back(vertex);
    edgesForVertex.insert(std::pair<V, std::vector<E>>(vertex, std::vector<E>()));
    verticesForVertex.insert(std::pair<V, std::vector<V>>(vertex, std::vector<V>()));

    return vertices.size() - 1;
  }

  /** \brief Adds an edge into this graph that binds two given vertices.
   *
   * <p>The method binds two vertices that belong to this graph with an edge. If any of the to
   * vertices is not found in the graph, it will be immediately inserted.
   * Since this graph is undirected, the order in which the two vertices
   * are given is unimportant</p>
   *
   * <p>This method forces the properties of a simple graph, i.e. it checks whether
   * the two given vertices are already connected. If so, a new edge is not added and
   * <strong>false</strong> value returned. Otherwise the method inserts the new edge and
   * returns <strong>true</strong>.
   *
   * @param firstVertex -index of the two vertexes to be connected
   * @param firstVertex -index of the two vertexes to be connected (the order is unimportant)
   * @param edge - the edge object
   * @return false if the two vertices are already connected by an edge; true otherwise
   * @see add_vertex(V)
   */
  bool add_edge(const core::index4 firstVertex, const core::index4  secondVertex, E &  edge) {

    if((firstVertex>=vertices.size())||(firstVertex>=vertices.size()))
      return false;
    V & v1 = vertices[firstVertex];
    V & v2 = vertices[secondVertex];
    if (!are_connected(v1, v2)) {
      verticesForEdge.insert(std::pair<E,std::pair<V,V>>(edge, std::pair<V,V>(v1, v2)));
      edges.push_back(edge);
      edgesForVertex.at(v1).push_back(edge);
      edgesForVertex.at(v2).push_back(edge);
      verticesForVertex.at(v1).push_back(v2);
      verticesForVertex.at(v2).push_back(v1);
      return true;
    }
    return false;
  }

  /** \brief Adds an edge into this graph that binds two given vertices.
   *
   * <p>The method binds two vertices that belong to this graph with an edge. If any of the to
   * vertices is not found in the graph, it will be immediately inserted.
   * Since this graph is undirected, the order in which the two vertices
   * are given is unimportant</p>
   *
   * <p>This method forces the properties of a simple graph, i.e. it checks whether
   * the two given vertices are already connected. If so, a new edge is not added and
   * <strong>false</strong> value returned. Otherwise the method inserts the new edge and
   * returns <strong>true</strong>.
   *
   * @param firstVertex - where an edge begins (the order is unimportant)
   * @param secondVertex - where an edge ends (the order is unimportant)
   * @param edge - the edge object
   * @return false if the two vertices are already connected by an edge; true otherwise
   */
  bool add_edge(V & firstVertex, V &  secondVertex, E &  edge) {

    if (!has_vertex(firstVertex)) add_vertex(firstVertex);

    if (!has_vertex(secondVertex)) add_vertex(secondVertex);

    if (!are_connected(firstVertex, secondVertex)) {
      verticesForEdge.insert(std::pair<E,std::pair<V,V>>(edge, std::pair<V,V>(firstVertex, secondVertex)));
      edges.push_back(edge);
      edgesForVertex.at(firstVertex).push_back(edge);
      edgesForVertex.at(secondVertex).push_back(edge);
      verticesForVertex.at(firstVertex).push_back(secondVertex);
      verticesForVertex.at(secondVertex).push_back(firstVertex);
      return true;
    }
    return false;
  }

  /** @brief Removes an edge between two vertices if it exists
   *
   * @param firstVertex - the first vertex
   * @param secondVertex - the second vertex
   * @return if the edge existed and was removed - returns true; false if there was no such edge
   */
  bool remove_edge(const V & firstVertex,const V &  secondVertex) {

    if(!are_connected(firstVertex,secondVertex)) return false;
    E e = get_edge(firstVertex, secondVertex);
    auto e_it = std::find(edges.begin(), edges.end(), e);
    edges.erase(e_it); // --- remove the edge from the vector of all edges of this graph
    verticesForEdge.erase(e); // --- remove the edge from edge-to-vertices map

    edgesForVertex.at(firstVertex).erase(std::find(edgesForVertex.at(firstVertex).begin(),edgesForVertex.at(firstVertex).end(),e));
    edgesForVertex.at(secondVertex).erase(std::find(edgesForVertex.at(secondVertex).begin(),edgesForVertex.at(secondVertex).end(),e));

    verticesForVertex.at(firstVertex).erase(std::find(verticesForVertex.at(firstVertex).begin(),verticesForVertex.at(firstVertex).end(),secondVertex));
    verticesForVertex.at(secondVertex).erase(std::find(verticesForVertex.at(secondVertex).begin(),verticesForVertex.at(secondVertex).end(),firstVertex));
    return true;
  }

  /** \brief returns the number of vertices in this graph.
   * @return the number of vertices
   */
  core::index4 count_vertices() const { return vertices.size(); }

  /** \brief Returns a being iterator for vertexes.
   * @return iterator that points prior the first vertex in this graph
   */
  typename std::vector<V>::iterator begin_vertex() { return vertices.begin(); }

  /**  \brief Returns an end iterator for vertexes.
   * @return iterator that points past the last vertex in this graph
   */
  typename std::vector<V>::iterator end_vertex() { return vertices.end(); }

  /** \brief Returns a being iterator for vertexes that are connected to <code>a_vertex</code>.
   * @return iterator that points prior the first vertex connected to <code>a_vertex</code>.
   */
  typename std::vector<V>::iterator begin_vertex(const V & a_vertex) { return verticesForVertex[a_vertex].begin(); }

  /** \brief Returns an end iterator for vertexes that are connected to <code>a_vertex</code>.
   * @return iterator that points past the last vertex connected to <code>a_vertex</code>.
   */
  typename std::vector<V>::iterator end_vertex(const V & a_vertex) { return verticesForVertex[a_vertex].end(); }

  /** \brief Returns a being iterator for vertexes.
 * @return iterator that points prior the first vertex in this graph
 */
  typename std::vector<V>::const_iterator cbegin_vertex() const { return vertices.cbegin(); }

  /**  \brief Returns an end iterator for vertexes.
   * @return iterator that points past the last vertex in this graph
   */
  typename std::vector<V>::const_iterator cend_vertex() const { return vertices.cend(); }

  /** \brief Returns a being iterator for vertexes that are connected to <code>a_vertex</code>.
   * @return iterator that points prior the first vertex connected to <code>a_vertex</code>.
   */
  typename std::vector<V>::const_iterator cbegin_vertex(const V & a_vertex) const { return verticesForVertex.at(a_vertex).cbegin(); }

  /** \brief Returns an end iterator for vertexes that are connected to <code>a_vertex</code>.
   * @return iterator that points past the last vertex connected to <code>a_vertex</code>.
   */
  typename std::vector<V>::const_iterator cend_vertex(const V & a_vertex) const { return verticesForVertex.at(a_vertex).cend(); }

  /** \brief Check the order of a given vertex.
   *
   * <p>The method counts how many edges points out from the given vertex. This is also the number of
   * vertices connected to this vertex. Note that loops also increase the order.</p>
   *
   * @param vertex
   *
   * @return the number of edges pointing out from the given vertex
   */
  core::index2 order(const V & vertex) const { return edgesForVertex.get(vertex).size(); }

  /**  \brief Returns the number of edges in this graph.
   * @return the number of edges
   */
  core::index4 countEdges() const { return edges.size(); }

  /**  \brief Returns a being iterator for edges.
   * @return iterator that points prior the first edge in this graph
   */
  typename std::vector<E>::iterator begin_edge() { return edges.begin(); }

  /**  \brief Returns an end iterator for edges.
   * @return iterator that points past the last edge in this graph
   */
  typename std::vector<E>::iterator end_edge() { return edges.end(); }

  /**  \brief Returns a being iterator for edges attached to a given vertex.
   * @return iterator that points prior the first edge of a given vertex in this graph
   */
  typename std::vector<E>::iterator begin_edge(const V & vertex) { return edgesForVertex[vertex].begin(); }

  /**  \brief Returns an end iterator for edges attached to a given vertex.
   * @return iterator that points past the last edge of a given vertex in this graph
   */
  typename std::vector<E>::iterator end_edge(const V & vertex) { return edgesForVertex[vertex].end(); }

  /**  \brief Returns a being const-iterator for edges.
   * @return iterator that points prior the first edge in this graph
   */
  typename std::vector<E>::const_iterator cbegin_edge() const { return edges.cbegin(); }

  /**  \brief Returns an end const-iterator for edges.
   * @return iterator that points past the last edge in this graph
   */
  typename std::vector<E>::const_iterator cend_edge() const { return edges.cend(); }

  /**  \brief Returns a being const-iterator for edges attached to a given vertex.
   * @return iterator that points prior the first edge of a given vertex in this graph
   */
  typename std::vector<E>::const_iterator cbegin_edge(const V & vertex) const { return edgesForVertex.at(vertex).cbegin(); }

  /**  \brief Returns an end const-iterator for edges attached to a given vertex.
   * @return iterator that points past the last edge of a given vertex in this graph
   */
  typename std::vector<E>::const_iterator cend_edge(const V & vertex) const { return edgesForVertex.at(vertex).cend(); }

  /**  \brief Returns an edge that connects the two given vertexes.
   */
  const E & get_edge(const V & vertex1,const V & vertex2) {
    for (const E & e : edgesForVertex[vertex1]) {
      const std::pair<V, V> & vv = verticesForEdge.at(e);
      if (((vv.first == vertex1) && (vv.second == vertex2)) || ((vv.first == vertex2) && (vv.second == vertex1)))
        return e;
    }
    throw std::runtime_error("Edge not found");
  }

  /**  \brief  Returns the two vertices that are connected by an edge.
   * @param edge - an edge that belongs to this graph
   * @return the TwoTuple object holding two vertices
   */
  std::pair<V, V> & get_edge_vertices(const E & edge) { return verticesForEdge.at(edge); }

  /** \brief Checks if two vertices are connected by at least one edge.
   *
   * @param vertex1 - a vertex from this graph
   * @param vertex2 - a partner to be checked
   * @return true if the two vertices are directly connected by a single edge
   */
  bool are_connected(const V &  vertex1,const V & vertex2) const {
    const std::vector<V> & v = verticesForVertex.at(vertex1);
    return (std::find(v.begin(), v.end(), vertex2) != v.end());
  }

protected:
  utils::Logger logger;
  std::vector<V> vertices;
  std::vector<E> edges;
  std::map<E, std::pair<V, V>> verticesForEdge;
  std::map<V, std::vector<E>> edgesForVertex;
  std::map<V, std::vector<V>> verticesForVertex;
};

}
}

#endif
