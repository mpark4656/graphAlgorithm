
#ifndef GRAPH_H
#define GRAPH_H

#include <utility>
#include <list>
#include <vector>
#include <iostream>
#include <limits>
#include <algorithm>

template <typename Type>
class Graph
{
private:
	/**
	 *	Data structure for the Node or vertex of the
	 *	graph.
	 *	Vertex class does not provide its own destructor
	 *	The garbage collection will be performed by the destructor of the Graph class.
	 *	Upon closer inspection, you will notice that garbage collection is not really
	 *	needed for the Vertex class because the memory addresses that the edges point to
	 *	are still used (they still exist in the graph) by the Graph class in the vertices
	 *	list. When the destructor of the Graph class is invoked, the memory will be freed.
	 */
	class Vertex
	{
	public:
		// Represents the state of this vertex.
		//	UNDISCOVERED - User never visited this vertex
		//	DISCOVERED - User visited this vertex, but hasn't yet visited all of its
		//			     adjacent vertices
		//	PROCESSED - User visited this vertex AND visited all of its adjacent vertices
		enum State {
			UNDISCOVERED ,
			DISCOVERED ,
			PROCESSED ,
		};
		
		/**
		 *	Constructor for Vertex
		 *	@param x Type data
		 */
		Vertex( Type x ): data{x}
		{}
		
		/**
		 *	Public Get function for Data
		 *	@return data Type data that this vertex contains
		 */
		Type getData() const
		{ return data; }
		
		/**
		 *	Public Get function for Degree
		 *	@return degree int outdegree of this vertex 
		 */
		int getDegree() const
		{ return edges.size(); }
		
		/**
		 *	Public Get function for weight of one of the edges in adjacency list
		 *	@param vertex Vertex* a pointer to the vertex
		 */
		double getWeight(Vertex* vertex) const
		{
			for(auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				if(itr->first == vertex) {
					return itr->second;
				}
			}
		}
		
		/**
		 *	This is used when traversing vertices in the graph
		 *	Specifies whether this vertex is Undiscovered.
		 *	@return true bool if this node or vertex has not been visited
		 */
		bool isUndiscovered() const
		{ return state == UNDISCOVERED; }
		
		/**
		 *	This is used when traversing vertices in the graph
		 *	Specifies whether this vertex has been discovered or not.
		 *	@return true bool if this node or vertex has been visited
		 */
		bool isDiscovered() const
		{ return state == DISCOVERED; }
		
		/**
		 *	This is used when traversing vertices in the graph
		 *	Specifies whether this vertex has been processed.
		 *	@return true bool if this vertex and all of its adjacent vertices have been
		 *		visited
		 */
		bool isProcessed() const
		{ return state == PROCESSED; }
		
		/**
		 *	Function to return the state of this vertex.
		 */
		State getState() const
		{ return state; }
		
		/**
		 *	Function to mark this vertex UNDISCOVERED.
		 */
		void markUndiscovered()
		{ state = UNDISCOVERED; }
		
		/**
		 *	Function to unmark this vertex DISCOVERED.
		 */
		void markDiscovered()
		{ state = DISCOVERED; }
		
		/**
		 *	Function to unmark this vertex PROCESSED.
		 */
		void markProcessed()
		{ state = PROCESSED; }
		
		/**
		 *	Public Set function for data
		 *	@param newData Type new value of data
		 */
		void setData( const Type &newData )
		{ data = newData; }
		
		/**
		 *	Public Set function for weight
		 *	@param vertex Vertex* defines which edge to update
		 *	@param newWeight double new value of weight
		 */
		void setWeight(Vertex* vertex , double newWeight)
		{
			// Find the iter at which the given vertex is found and set the weight
			// (second) to the new weight.
			for(auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				if(itr->first == vertex) {
					itr->second = newWeight;
				}
			}
		}
		
		/**
		 *	Add an edge for this vertex
		 *	This function assumes that order in the adjacency list does not matter.
		 *	(And usually, order does not matter in the adjacency list)
		 *	If the graph is simple and duplicate edges must be avoided, it is the
		 *	responsibility of the Graph class to verify that there is no duplicate
		 *	before adding a new edge.
		 *	@param toAdd Vertex* a pointer to the destination vertex to add
		 *	@param weight double weight of the edge
		 */
		void addEdge(Vertex* toAdd , double weight)
		{
			edges.push_back( std::pair<Vertex*,double>{toAdd , weight} );
		}
		
		/**
		 *	Modify an existing edge for this vertex.
		 *	This function assumes that order in the adjacency list does not matter.
		 *	(And usually, order does not matter in the adjacency list)	 
		 *	@param toModify Vertex* a pointer to the destination vertex to modify
		 *	@param weight double weight of the existing edge
		 *	@param newWeight double new value of the weight for this edge
		 */
		void modifyEdge(Vertex* toModify , double weight , double newWeight)
		{
			// Check to see if this vertex already exists in the adjacency list
			// If it already exists, remove it from the list and re-add it
			// (The new vertex may have updated information)
			// Notice that this compares the pointer variables or the memory
			// "address" of the vertices.
			for(auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				if(itr->first == toModify && itr->second == weight) {
					edges.erase(itr);
					edges.push_back( std::pair<Vertex*,double>{toModify , newWeight} );
					break;
				}
			}
		}
		
		/**
		 *	Erase an edge from this vertex.
		 *	@param toErase Vertex* a pointer to the vertex to erase
		 */
		void eraseEdge(Vertex* toErase , double weight)
		{
			for(auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				if(itr->first == toErase && itr->second == weight) {
					edges.erase(itr);
					break;
				}
			}
		}
		
		/**
		 *	Clear all edges for this vertex.
		 */
		void clearAllEdges()
		{
			edges.clear();
		}
		
		/**
		 *	Public function to return a vector(array) of data
		 *	in this vertex's adjacency list
		 *	@return result std::vector<Type> list of data in the adjacency list
		 */
		std::vector<Type> dataAdjacencyList() const
		{
			std::vector<Type> result;
			
			for( auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				result.push_back(itr->first->data);
			}
			
			return result;
		}
		
		/**
		 *	Public function to return a vector(array) of weights of edges
		 *	in this vertex's adjacency list
		 *	@return result std::vector<double> list of weights in the adjacency list
		 */
		std::vector<double> weightAdjacencyList() const
		{
			std::vector<double> result;
			
			for( auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				result.push_back(itr->second);
			}
			
			return result;
		}
		
		/**
		 *	Finds out whether the given vertex is in the adjacency list for this vertex
		 *	@param toFind Vertex* vertex to find
		 *	@return true bool if the vertex is in the adjacency list of this vertex
		 */
		bool isInAdjacencyList(Vertex* toFind) const
		{
			for( auto itr = edges.begin() ; itr != edges.end() ; itr++) {
				if(itr->first == toFind) {
					return true;
				}
			}
			
			return false;
		}
		
	private:
		// This is a unique key in the graph. There cannot be vertices with duplicate keys
		// The key provide its own operator==
		Type data;
		
		// Used when traversing the graph.
		State state = UNDISCOVERED;
		
		// List of Pairs of next Node and
		// corresponding weight
		std::list<std::pair<Vertex* , double>> edges;
		
	}; // End of Class Vertex

public:
	// This is the same as State enum in the Vertex class
	// However, these enum variables are available for the user to use.
	// The Vertex class is private and it's hidden from the user
	enum VertexState {
		NON_EXISTENT ,
		UNDISCOVERED ,
		DISCOVERED ,
		PROCESSED ,
	};
	
	/**
	 *	Default Constructor
	 */
	Graph();

	/**
	 *	Constructor that accepts 2 parameters
	 *	@param directed bool Is the graph directed?
	 */
	Graph( bool isDirected );
	
	/**
	 *	Copy Constructor
	 */
	Graph( const Graph &toCopy );
	
	/**
	 *	Move Constructor
	 */
	Graph( Graph &&toCopy );
	
	/**
	 *	Assignment Operator
	 */
	Graph & operator= ( const Graph &rhs );
	
	/**
	 *	Move Assignment Operator
	 */
	Graph & operator= ( Graph &&rhs );
	
	/**
	 *	Destructor
	 */
	~Graph();
	
	/**
	 *	Add a vertex to this graph
	 *	@param x Type the vertex to add
	 */
	void addVertex( const Type &x );

	/**
	 *	Remove the specified vertex.
	 *	This will also remove ALL adjacencies that involve this vertex. Nothing will have
	 *	relationship with this vertex.
	 *	@param x Type the vertex to remove
	 */
	void removeVertex( const Type &x );
	
	/**
	 *	Change the datum at the vertex
	 *	@param x Type the vertex to change
	 *	@param y Type the new vertex
	 */
	void modifyVertex( const Type &x , const Type &newData );
	
	/**
	 *	Add an edge. How the edge is added really depends on whether
	 *	the graph is directed or undirected. For a directed graph, the adjacency is
	 *	added once. For undirected graph, the adjacency is added twice.
	 *	If the x and y don't exist, new vertices will be created.
	 *	This function is called without the weight parameter (Weight will be 0 by default)
	 *	User will call this function if they aren't really concerned with the weight.
	 *	@param x Type Represents "From" vertex
	 *	@param x Type Represents "To" vertex
	 */
	void addEdge( const Type &x , const Type &y );

	/**
	 *	Add an edge. How the edge is added really depends on whether
	 *	the graph is directed or undirected. For a directed graph, the adjacency is
	 *	added once. For undirected graph, the adjacency is added twice.
	 *	If the x and y don't exist, new vertices will be created.
	 *	Calls internal (private) functions to handle the actual addition of the edge.
	 *	@param x Type Represents "From" vertex
	 *	@param x Type Represents "To" vertex
	 */
	void addEdge( const Type &x , const Type &y , double weight );
	
	/**
	 *	Modify the weight of currently existing edge.
	 *	Actually, calling addEdge() again with the same parameters will
	 *	overwrite the existing edge with whatever new weight value is given and
	 *	this function does exactly that (Calling addEdge() again).
	 *	If the given data do not exist, nothing will happen.
	 *	@param x Type the existing "From" vertex
	 *	@param y Type the existing "To" vertex
	 *	@param weight double new weight value of this edge or relationship
	 */
	void modifyEdge( const Type &x , const Type &y , double weight , double newWeight );
	
	/**
	 *	Remove the specified edge from this graph. This will NOT remove the vertices
	 *	themselves. It just removes the relationship. If you need to remove the vertices,
	 *	call removeVertex() function.
	 *	User can call this function if they aren't really concerned with the weight of the
	 *	edge. It will be assumed that the weight is simply 0 (If you called addEdge 
	 *	without the weight parameter, all edges will have 0 weight).
	 *	@param x Type the "From" vertex
	 *	@param y Type the "To" vertex
	 */
	void removeEdge( const Type &x , const Type &y );

	/**
	 *	Remove the specified edge from this graph. This will NOT remove the vertices
	 *	themselves. It just removes the relationship. If you need to remove the vertices,
	 *	call removeVertex() function.
	 *	@param x Type the "From" vertex
	 *	@param y Type the "To" vertex
	 *	@param weight double the weight of the edge
	 */
	void removeEdge( const Type &x , const Type &y , double weight );
	
	void removeOutDegreeOf( const Type &x );
	
	void removeInDegreeOf( const Type &x );
	
	void isolate( const Type &x );
	
	/**
	 *	Return the out-degree of the specified vertex.
	 *	@param x Type vertex
	 *	@return degree int degree of the vertex. -1 if the vertex doesn't exist
	 */
	int outDegreeOf( const Type &x ) const;
	
	/**
	 *	Return the in-degree of the specified vertex.
	 *	For undirected graph, in-degree == out-degree
	 *	@param x Type vertex
	 *	@return degree int degree of the vertex. -1 if the vertex doesn't exist
	 */
	int inDegreeOf( const Type &x ) const;
	
	/**
	 *	Return the weight of the specified edge
	 *	If there are multiple edges with same vertices, then the weight of the first edge
	 *	on the adjacency list will be returned.
	 *	Return -std::numeric_limit<double>::max if the edge doesn't exist.
	 */
	double weightOf( const Type &x , const Type &y ) const;
	
	/**
	 *	Return the total number of vertices in this graph.
	 */
	int totalVertices() const;
	
	/**
	 *	Return the total number of edges in this graph.
	 *	The calculation method is a little different depending on
	 *	whether the graph is directed or undirected.
	 */
	int totalEdges() const;
	
	/**
	 *	Return the list of all vertices in the graph
	 *	The returned data structure is std::vector
	 *	@return result std::vector<Type> list of all vertices in the graph
	 */
	std::vector<Type> verticesList() const;
	
	/**
	 *	Function to return the adjacency list of the given vertex
	 *	The returned data structure is std::vector
	 *	@param x Type vertex
	 *	@return adjacencyList std::vector<Type> list of adjacent vertices
	 */
	std::vector<Type> adjacentVertexList( const Type &x ) const;
	
	/**
	 *	Function to return the weight of adjacent vertices of the given vertex
	 *	The returned data structure is std::vector
	 *	@param x Type vertex
	 *	@return adjacentWeights std::vector<double> weights of adjacent vertices
	 */
	std::vector<double> adjacentWeightList( const Type &x ) const;
	
	/**
	 *	Function to mark a vertex as Undiscovered
	 *	@param x Type Vertex
	 */
	void undiscover( const Type &x );
	
	/**
	 *	Function to mark a vertex as Discovered
	 *	@param x Type Vertex
	 */
	void discover( const Type &x );
	
	/**
	 *	Function to mark a vertex as Processed
	 *	@param x Type Vertex
	 */
	void process( const Type &x );
	
	/**
	 *	Function to return true if this vertex is UNDISCOVERED
	 *	@param x Type vertex
	 */
	bool isUndiscovered( const Type &x ) const;
	
	/**
	 *	Function to return true if this vertex is DISCOVERED
	 *	@param x Type vertex
	 */
	bool isDiscovered( const Type &x ) const;
	
	/**
	 *	Function to return true if this vertex is PROCESSED
	 *	@param x Type vertex
	 */	
	bool isProcessed( const Type &x ) const;
	
	void undiscoverAll();
	
	void discoverAll();
	
	void processAll();
	
	void clearSelfAdjacency( const Type &x );
	
	void removeAllVertices();
	
	void removeAllEdges();
	
	/**
	 *	Returns the state of the given vertex.
	 *	@param x Type vertex
	 *	@return state VertexState state of the given vertex (enum)
	 */
	VertexState getState( const Type &x ) const;
	
	/**
	 *	The function to test whether the given vertex exists.
	 *	@param x Type vertex
	 *	@return true bool if the vertex is in the graph
	 */
	bool exists( const Type &x ) const;
	
	/**
	 *	The function to test whether the given edge exists.
	 *	@param x Type the "From" vertex
	 *	@param y Type the "To" vertex
	 *	@return true bool if the edge exists
	 */
	bool exists( const Type &x , const Type &y ) const;
	
	/**
	 *	Returns true if the graph is directed. Returns false otherwise
	 */
	bool isDirected() const;
	
	/**
	 *	Display the vertices and the corresponding adjacency lists
	 */
	void display( std::ostream &out ) const;
	
private:
	// List of Vertices in this graph
	std::list<Vertex*> vertices;
	
	// Specifies whether this graph is directed or not
	bool directed;
	
	/**
	 *	Returns pointer variable for the given data
	 *	Returns nullptr if the given data does not yet exist in this graph.
	 *	@param x Type data value
	 *	@return ptr Vertex* pointer to the vertex that holds the data
	 */
	Vertex* getPointer( const Type &x ) const;
	
}; // End of Class Graph


/**
 *	Returns pointer variable for the given data
 *	Returns nullptr if the given data does not yet exist in this graph.
 *	@param x Type data value
 *	@return ptr Vertex* pointer to the vertex that holds the data
 */
template <typename Type>
typename Graph<Type>::Vertex* Graph<Type>::getPointer( const Type &x ) const
{
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		if( (*itr)->getData() == x ) {
			return *itr;
		}
	}
	
	return nullptr;
}

/**
 *	Default Constructor
 */
template <typename Type>
Graph<Type>::Graph()
: directed{false}
{}
	
/**
 *	Constructor that accepts a parameter
 *	@param directed bool Is the graph directed?
 */
template <typename Type>
Graph<Type>::Graph( bool isDirected )
: directed{isDirected}
{}

/**
 *	Copy Constructor
 */
template <typename Type>
Graph<Type>::Graph( const Graph &toCopy )
{
	// During the copying process, we want the graph to add only one edge
	// So, the directed flag is temporarily set to true
	this->directed = true;
	
	std::vector<Type> vertices = toCopy.verticesList();
	
	for( int i = 0 ; i < vertices.size() ; i++ ) {
		std::vector<Type> adjacentVertices = 
			toCopy.adjacentVertexList( vertices[i] );
			
		if( adjacentVertices.size() == 0 ) {
			this->addVertex( vertices[i] );
		}
		else {
			std::vector<double> adjacentWeights =
				toCopy.adjacentWeightList( vertices[i] );
			
			for( int j = 0 ; j < adjacentVertices.size() ; j++ ) {
				this->addEdge( vertices[i] , adjacentVertices[j] , adjacentWeights[j] );
			}
		}
	}
	
	// Specifiy whether this graph should be directed or undirected
	this->directed = toCopy.isDirected();
}

/**
 *	Move Constructor
 */
template <typename Type>
Graph<Type>::Graph( Graph &&toCopy )
: directed{ toCopy.directed }
{
	toCopy.vertices.swap( this->vertices );
}

/**
 *	Assignment Operator
 */
template <typename Type>
Graph<Type> & Graph<Type>::operator= ( const Graph &rhs )
{
	Graph<Type> copy = rhs;
	std::swap( *this , copy );
	//std::swap( this->directed , copy.directed );
	//std::swap( this->vertices , copy.vertices );
	return *this;
}

/**
 *	Move Assignment Operator
 */
template <typename Type>
Graph<Type> & Graph<Type>::operator= ( Graph &&rhs )
{
	std::swap( this->directed , rhs.directed );
	std::swap( this->vertices , rhs.vertices );
	return *this;
}

/**
 *	Destructor
 */
template <typename Type>
Graph<Type>::~Graph()
{
	typename std::list<Graph<Type>::Vertex*>::iterator itr = vertices.begin();
	
	while( itr != vertices.end() ) {
		auto toDelete = itr;
		itr++;
		delete *toDelete;
	}
}

/**
 *	Add a vertex to this graph
 *	@param x Type the vertex to add
 */
template <typename Type>
void Graph<Type>::addVertex( const Type &x )
{
	// Check if this vertex already exists
	if( exists(x) ) {
		return;
	}
	
	vertices.push_back( new Vertex{x} );
}

/**
 *	Remove the specified vertex.
 *	@param x Type vertex
 */
template <typename Type>
void Graph<Type>::removeVertex( const Type &x )
{
	// Check if the vertex is not in the graph, take no action.
	if( !exists(x) ) {
		return;
	}
	
	// Store the pointer in a temporary variable
	Vertex* toDelete = getPointer(x);
	
	// First, remove this vertex from the adjacency list of all vertices
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		(*itr)->eraseEdge(toDelete , (*itr)->getWeight(toDelete) );
	}
	
	// Second, remove the vertex from the list of vertices (From this graph)
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		if( (*itr)->getData() == x ) {
			vertices.erase(itr);
			break;
		}
	}
	
	// Deallocate the memory
	delete toDelete;
}

/**
 *	Change the datum at the vertex
 *	@param x Type the existing vertex
 *	@param y Type the new vertex datum
 */
template <typename Type>
void Graph<Type>::modifyVertex( const Type &x , const Type &newData )
{
	// Check if the vertex is not in the graph, take no action.
	if( !exists(x) ) {
		return;
	}
	
	// Check to see if x == newData. If so, then user may be trying
	// to change underlying information.
	if( x == newData ) {
		; // Allow the program to execute the next statement after if-else
	}
	// If x and newData have different identifiers (keys), make sure there is no other
	// vertex with the same key as newData, then update x
	else {
		// If there is a vertex with the same identifier as newData, then
		// we can't update x with newData since that will results in duplicate keys
		// in the graph. Take no further action.
		for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
			if( (*itr)->getData() == newData ) {
				return;
			}
		}
	}
	
	// If the control has come this far, it means the newData has a unique key and it
	// can be added to the graph OR the key for x and newData is the same. 
	// Go ahead and change x to the newData.
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		if( (*itr)->getData() == x ) {
			(*itr)->setData(newData);
			return;
		}
	}
}

/**
 *	Add an edge. How the edge is added really depends on whether
 *	the graph is directed or undirected. For a directed graph, the adjacency is
 *	added once. For undirected graph, the adjacency is added twice.
 *	If the x and y don't exist, new vertices will be created.
 *	This function is called without the weight parameter(Weight will be 0 by default).
 *	User can call this function if they aren't really concerned with the weight of the
 *	edge.
 *	@param x Type Represents "From" vertex
 *	@param x Type Represents "To" vertex
 */
template <typename Type>
void Graph<Type>::addEdge( const Type &x , const Type &y )
{
	// Add edge with 0 weight.
	addEdge(x , y , 0);
}

/**
 *	Add an edge. How the edge is added really depends on whether
 *	the graph is directed or undirected. For a directed graph, the adjacency is
 *	added once. For undirected graph, the adjacency is added twice.
 *	If the x and y don't exist, new vertices will be created.
 *	Calls internal (private) functions to handle the actual addition of the edge.
 *	@param x Type Represents "From" vertex
 *	@param x Type Represents "To" vertex
 */
template <typename Type>
void Graph<Type>::addEdge( const Type &x , const Type &y , double weight )
{
	// If the edge(x , y) with the same weight already exists, do nothing.
	if( exists(x , y) && weightOf(x , y) == weight ) {
		return;
	}
	
	// If x doesn't exist, create a new vertex.
	if( !exists(x) ) {
		addVertex(x);
	}
	
	// If y doesn't exist, create a new vertex.
	if( !exists(y) ) {
		addVertex(y);
	}
	
	// Add the edge
	if(isDirected()) {
		getPointer(x)->addEdge( getPointer(y) , weight );
	}
	else {	// Undirected Graph
		getPointer(x)->addEdge( getPointer(y) , weight );
		getPointer(y)->addEdge( getPointer(x) , weight );
	}
}
	
/**
 *	Modify the weight of currently existing edge.
 *	Actually, calling addEdge() again with the same parameters will
 *	overwrite the existing edge with whatever new weight value is given and
 *	this function does exactly that (Calling addEdge() again).
 *	If the given data do not exist, nothing will happen.
 *	@param x Type the existing "From" vertex
 *	@param y Type the existing "To" vertex
 *	@param weight double weight value of the existing edge
 *	@param newWeight double new value of the weight
 */
template <typename Type>
void Graph<Type>::modifyEdge( 
	const Type &x ,
	const Type &y ,
	double weight ,
	double newWeight
)
{
	// If vertex x or y doesn't exist, take no further action
	if( !exists(x) || !exists(y) ) {
		return;
	}
	
	// Modify the edge
	if(isDirected()) {
		getPointer(x)->modifyEdge(getPointer(y) , weight , newWeight);
	}
	else {	// Undirected Graph
		getPointer(x)->modifyEdge(getPointer(y) , weight , newWeight);
		getPointer(y)->modifyEdge(getPointer(x) , weight , newWeight);
	}
}

/**
 *	Remove the specified edge from this graph. This will NOT remove the vertices
 *	themselves. It just removes the relationship. If you need to remove the vertices,
 *	call removeVertex() function.
 *	User can call this function if they aren't really concerned with the weight of the
 *	edge. It will be assumed that the weight is simply 0 (If you called addEdge without
 *	the weight parameter, all edges will have 0 weight).
 *	@param x Type the "From" vertex
 *	@param y Type the "To" vertex
 */
template <typename Type>
void Graph<Type>::removeEdge( const Type &x , const Type &y )
{
	removeEdge( x , y , 0 );
}

/**
 *	Remove the specified edge from this graph. This will NOT remove the vertices
 *	themselves. It just removes the relationship. If you need to remove the vertices,
 *	call removeVertex() function.
 *	@param x Type the "From" vertex
 *	@param y Type the "To" vertex
 *	@param weight double the weight of the edge
 */
template <typename Type>
void Graph<Type>::removeEdge( const Type &x , const Type &y , double weight )
{
	// If the given edge doens't exist, take no further action
	if( !exists(x , y) ) {
		return;
	}
	
	if(isDirected()) {
		getPointer(x)->eraseEdge( getPointer(y) , weight );
	}
	else {	// Undirected
		getPointer(x)->eraseEdge( getPointer(y) , weight );
		getPointer(y)->eraseEdge( getPointer(x) , weight );
	}
}

/**
 *	Return the out-degree of the vertex.
 *	@param x Type vertex
 *	@return degree int degree of the vertex. -1 if the vertex doesn't exist
 */
template <typename Type>
int Graph<Type>::outDegreeOf( const Type &x ) const
{
	// If the vertex doesn't exist, return -1
	if( !exists(x) ) {
		return -1;
	}
	
	return getPointer(x)->getDegree();
}

/**
 *	Return the in-degree of the given vertex.
 *	For undirected graph, in-degree == out-degree
 *	@param x Type vertex
 *	@return degree int degree of the vertex. -1 if the vertex doesn't exist
 */
template <typename Type>
int Graph<Type>::inDegreeOf( const Type &x ) const
{
	// If the vertex doesn't exist, return -1
	if( !exists(x) ) {
		return -1;
	}
	
	int degree = 0;
	
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		if( (*itr)->isInAdjacencyList(getPointer(x)) ) {
			degree++;
		}
	}
	
	return degree;
}

/**
 *	Return the weight of the specified edge
 *	If there are multiple edges with same vertices, then the weight of the first edge
 *	on the adjacency list will be returned.
 *	Return -std::numeric_limits<double>::max() if the edge doesn't exist.
 */
template <typename Type>
double Graph<Type>::weightOf( const Type &x , const Type &y ) const
{
	if( !exists(x , y) ) {
		return -std::numeric_limits<double>::max();
	}
	
	return getPointer(x)->getWeight( getPointer(y) );
}

/**
 *	Return the total number of vertices in this graph.
 */
template <typename Type>
int Graph<Type>::totalVertices() const
{
	return vertices.size();
}

/**
 *	Return the total number of edges in this graph.
 *	The calculation method is a little different depending on
 *	whether the graph is directed or undirected.
 */
template <typename Type>
int Graph<Type>::totalEdges() const
{
	int counts = 0;
	
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		counts += (*itr)->getDegree();
	}
		
	if( !isDirected() ) {
		counts /= 2;
	}
	
	return counts;
}

/**
 *	Return the list of all vertices in the graph
 *	The returned data structure is std::vector
 *	@return result std::vector<Type> list of all vertices in the graph
 */
template <typename Type>
std::vector<Type> Graph<Type>::verticesList() const
{
	std::vector<Type> result;
	
	for( auto itr = vertices.begin() ; itr != vertices.end() ; itr++ ) {
		result.push_back( (*itr)->getData() );
	}
	
	return result;
}
	
/**
 *	Function to return the adjacency list of the given vertex
 *	The returned data structure is std::vector
 *	@param x Type vertex
 *	@return adjacencyList std::vector<Type> list of adjacent vertices
 */
template <typename Type>
std::vector<Type> Graph<Type>::adjacentVertexList( const Type &x ) const
{
	if( exists(x) ) {
		return getPointer(x)->dataAdjacencyList();
	}
	else {
		//Return an empty vector if the provided vertex does not exist
		return std::vector<Type>{};
	}
}

/**
 *	Function to return the weight of adjacent vertices of the given vertex
 *	The returned data structure is std::vector
 *	@param x Type vertex
 *	@return adjacentWeights std::vector<double> weights of adjacent vertices
 */
template <typename Type>
std::vector<double> Graph<Type>::adjacentWeightList( const Type &x ) const
{
	if( exists(x) ) {
		return getPointer(x)->weightAdjacencyList();
	}
	else {
		//Return an empty vector if the provided vertex does not exist
		return std::vector<double>{};
	}
	
}

/**
 *	Function to mark a vertex as Undiscovered
 *	@param x Type Vertex
 */
template <typename Type>
void Graph<Type>::undiscover( const Type &x )
{
	if( exists(x) ) {
		getPointer(x).markUndiscovered();
	}
}

/**
 *	Function to mark a vertex as Discovered
 *	@param x Type Vertex
 */
template <typename Type>
void Graph<Type>::discover( const Type &x )
{
	if( exists(x) ) {
		getPointer(x).markDiscovered();
	}
}

/**
 *	Function to mark a vertex as Processed
 *	@param x Type Vertex
 */
template <typename Type>
void Graph<Type>::process( const Type &x )
{
	if( exists(x) ) {
		getPointer(x).markProcessed();
	}
}

/**
 *	Function to return true if this vertex is UNDISCOVERED
 *	@param x Type vertex
 */
template <typename Type>
bool Graph<Type>::isUndiscovered( const Type &x ) const
{
	if( exists(x) ) {
		return getPointer(x)->isUndiscovered();
	}
	else {
		// SHOULD REALLY THROW EXCEPTION HERE
		return false;
	}
}

/**
 *	Function to return true if this vertex is DISCOVERED
 *	@param x Type vertex
 */
template <typename Type>
bool Graph<Type>::isDiscovered( const Type &x ) const
{
	if( exists(x) ) {
		return getPointer(x)->isDiscovered();
	}
	else {
		// SHOULD REALLY THROW EXCEPTION HERE
		return false;
	}
}

/**
 *	Function to return true if this vertex is PROCESSED
 *	@param x Type vertex
 */
template <typename Type>
bool Graph<Type>::isProcessed( const Type &x ) const
{
	if( exists(x) ) {
		return getPointer(x)->isProcessed();
	}
	else {
		// SHOULD REALLY THROW EXCEPTION HERE
		return false;
	}
}

template <typename Type>
void Graph<Type>::undiscoverAll()
{
	for( int i = 0 ; i < vertices.size() ; i++ ) {
		vertices[i]->markUndiscovered();
	}
}

/**
 *	Returns the state of the given vertex.
 *	@param x Type Vertex
 *	@return state VertexState state of the given vertex (enum)
 */
template <typename Type>
typename Graph<Type>::VertexState Graph<Type>::getState( const Type &x ) const
{
	// If the provided vertex doesn't exist, return NON_EXISTENT
	if( !exists(x) ) {
		return Graph<Type>::NON_EXISTENT;
	}
	
	if( getPointer(x)->getState() == Vertex::UNDISCOVERED ) {
		return Graph<Type>::UNDISCOVERED;
	}
	else if( getPointer(x)->getState() == Vertex::DISCOVERED ) {
		return Graph<Type>::DISCOVERED;
	}
	else {
		return Graph<Type>::PROCESSED;
	}
}

/**
 *	The function to test whether the given vertex exists
 *	@param x Type the Vertex to test
 *	@return true bool if the vertex is in the graph
 */
template <typename Type>
bool Graph<Type>::exists( const Type &x ) const
{
	if( getPointer(x) == nullptr ) {
		return false;
	}
	
	return true;
}

/**
 *	The function to test whether the given edge exists.
 *	@param x Type the "From" vertex
 *	@param y Type the "To" vertex
 *	@return true bool if the edge exists
 */
template <typename Type>
bool Graph<Type>::exists( const Type &x , const Type &y ) const
{
	// If vertex x or vertex y doesn't exist, then obviously, the edge doesn't exist
	if( !exists(x) || !exists(y) ) {
		return false;
	}
	
	return getPointer(x)->isInAdjacencyList( getPointer(y) );
}

/**
 *	Returns true if the graph is directed. Returns false otherwise
 */
template <typename Type>
bool Graph<Type>::isDirected() const
{
	return directed;
}

/**
 *	Display the vertices and the corresponding adjacency lists
 */
template <typename Type>
void Graph<Type>::display( std::ostream &out ) const
{
	// Specify if this is a directed graph
	out << "Directed = ";
	
	if(isDirected()) {
		out << "true\n";
	}
	else {
		out << "false\n";
	}
	
	// List the vertices and their adjacency lists
	int index = 1;
	
	for(auto itr = vertices.begin() ; itr != vertices.end() ; itr++) {
		// Display the index
		out << index << ". ";
		index++;
		
		out << (*itr)->getData() << "\n";
		
		std::vector<Type> data = (*itr)->dataAdjacencyList();
		std::vector<double> weights = (*itr)->weightAdjacencyList();
		
		for( int i = 0 ; i < data.size() ; i++ ) {
			out << " " << data[i] << " (" << weights[i] << ") \n";
		}
	}
}

#endif



// End of File
