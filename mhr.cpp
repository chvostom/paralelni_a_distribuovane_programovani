#include <iostream>
#include <string.h>
#include <cmath>
#include <unordered_map>
#include <set>
#include <string>
#include <iomanip>
#include <vector>
#include <deque>
#include <chrono>

using namespace std;
using namespace chrono;

/*==========================================CONSTANTS===========================================*/
/*
 * MAX_THREADS - count of threads 
 * ENOUGH_STATES - potřebný počet stavů, které musí mít STL kontejner přidělený vláknům
 */

#define MAX_THREADS 8
#define ENOUGH_STATES 800
/*================================GLOBAL VARIABLE - READ & WRITE================================*/
/*
 * minimalPrice - global variable to store the price of minimal edge cut
 * minimalCutState - global vector which describes both sets with vertices
 */

double minimalPrice = 0.0;
vector<bool> minimalCutState;

/*=========================================RESULT class=========================================*/
/*
 * Class which represents result of one dataset and contains:
 * price - price of minimal edge cut
 * complexity - number of function solve calls 
 * time - time of whole solve process in seconds
 * constructor without parameters
 * constructor with parameters
 * getters
 */

class Result{
public:
	Result( );
	Result( double price, double complexity, double time );
	double getPrice( );
	double getComplexity( );
	double getTime( );
private:
	double price;
	double complexity;
	double time;
};

/*----------------------------------------------------------------------------------------------*/

Result::Result( ){
	this->price = 0.0;
	this->complexity = 0.0;
	this->time = 0.0;
}

/*----------------------------------------------------------------------------------------------*/

Result::Result( double price, double complexity, double time ){
	this->price = price;
	this->complexity = complexity;
	this->time = time;
}

/*----------------------------------------------------------------------------------------------*/

double Result::getPrice( ) {
	return this->price;
}

/*----------------------------------------------------------------------------------------------*/

double Result::getComplexity( ) {
	return this->complexity;
}

/*----------------------------------------------------------------------------------------------*/

double Result::getTime( ) {
	return this->time;
}

/*==========================================EDGE class==========================================*/
/*
 * Class which represents edge between two vertices and contains:
 * startNode - index of first node
 * endNode - index of second node
 * price - price value of the edge
 * constructor with parameters
 * getters
 */

class Edge{
public:
	Edge( int fromNode, int toNode, double price );
	int getStartNode( );
	int getEndNode( );
	double getPrice( );
private:
	int startNode;
	int endNode;
	double price;
};

/*----------------------------------------------------------------------------------------------*/

Edge::Edge( int fromNode, int toNode, double price ) {
	this->startNode = fromNode;
	this->endNode = toNode;
	this->price = price;
}

/*----------------------------------------------------------------------------------------------*/

int Edge::getStartNode( ) {
	return this->startNode;
}

/*----------------------------------------------------------------------------------------------*/
	
int Edge::getEndNode( ) {
	return this->endNode;
}

/*----------------------------------------------------------------------------------------------*/

double Edge::getPrice( ) {
	return this->price;
}

/*=========================================GRAPH class==========================================*/
/*
 * Class that represents whole graph and contains:
 * nodeCount - count of vertices
 * edgeCount - count of edges
 * graph - vector which contains vector with edges for all vertices
 * constructor with parameters
 * addEdge - void which add one edge to the graph
 * getEdgesToNode - function which return vector with edges to concrete vertex of graph
 * loadGraph - void to load input
 */

class Graph{
public:
	Graph( int nodes, int edges );
	void addEdge( Edge edge );
	vector<Edge> getEdgesToNode( int node );
	void loadGraph( );
	vector<vector<Edge>> graph;
private:
	int nodeCount;
	int edgeCount;
};

/*----------------------------------------------------------------------------------------------*/

Graph::Graph( int nodes, int edges ) {
	for ( int i = 0 ; i < nodes ; i++ ){
		vector<Edge> vectorTmp;
		this->graph.push_back( vectorTmp );
	}
	this->nodeCount = nodes;
	this->edgeCount = edges;
	this->loadGraph( );
}

/*----------------------------------------------------------------------------------------------*/

void Graph::addEdge( Edge edge ){
	this->graph[edge.getEndNode( )].push_back( edge );
}

/*----------------------------------------------------------------------------------------------*/

vector<Edge> Graph::getEdgesToNode( int node ) {
	return this->graph[node];
}

/*----------------------------------------------------------------------------------------------*/

void Graph::loadGraph( ) {
	int startNode;
	int endNode;
	double price;
	for ( int i = 0 ; i < nodeCount * edgeCount / 2 ; i++ ) {
		cin >> startNode;
		cin >> endNode;
		cin >> price;
		minimalPrice += price;
		this->graph[endNode].push_back( Edge( startNode, endNode, price ) );
	}
	minimalPrice += 1.0;
}

/*=================================GLOBAL VARIABLE - READ ONLY==================================*/
/* 
 * graph - instance of Graph class which has information about graph (vertices and edges)
 * nodeCount - count of vertices
 * edgeCount - count of edges
 * conditionNumber - count of vertices that the set X must contain
 */

Graph graph = Graph( 0, 0 );
int nodeCount = 0;
int edgeCount = 0;
int conditionNumber = 0;

/*=========================================STATE class==========================================*/
/*
 * Class that represents concrete state of BFS algorithm and contains:
 * actualNode - actual vertex of graph 
 * actualPrice - actual price of subcut
 * actualX - count of vertices in set X
 * state - boolean vector that has information about content of sets X and Y
 * constructor with and without parameters
 * getNextBFSStates 
 	- function that return deque with descendants 
	- if state has no descendants then return deque with self
 */

class BFSState{
public:
	BFSState( );
	BFSState( int actualNode, double actualPrice, int actualX, vector<bool> state );
	deque<BFSState> getNextBFSStates( );
	int actualNode; 
	double actualPrice; 
	int actualX; 
	vector<bool> state;
};

/*----------------------------------------------------------------------------------------------*/

BFSState::BFSState( ){
	this->actualNode = 0;
	this->actualPrice = 0.0;
	this->actualX = 0;
}

/*----------------------------------------------------------------------------------------------*/

BFSState::BFSState( int actualNode, double actualPrice, int actualX, vector<bool> state ){
	this->actualNode = actualNode;
	this->actualPrice = actualPrice;
	this->actualX = actualX;
	this->state = state;
}

/*----------------------------------------------------------------------------------------------*/

deque<BFSState> BFSState::getNextBFSStates( ) {
	deque<BFSState> nextStates;
	if ( this->actualNode == ( nodeCount - 1 ) ) {
		nextStates.push_back( *this );
	}
	else{
		state.push_back( false );
		int actualY = actualNode - actualX;
		double actualPriceX = actualPrice;
		double actualPriceY = actualPrice;
		for ( int i = 0 ; i < ( int ) graph.graph[actualNode].size( ) ; i++ ){
			if( state[graph.graph[actualNode][i].getStartNode( )] ){
				actualPriceY += graph.graph[actualNode][i].getPrice( );
			}
			else {
				actualPriceX += graph.graph[actualNode][i].getPrice( );
			}
		}
		state[actualNode] = true;
		actualX++;
		if( ( nodeCount - actualY >= conditionNumber ) && ( actualX <= conditionNumber ) && ( actualPriceX < minimalPrice ) ) {
			nextStates.push_back( BFSState( actualNode + 1, actualPriceX, actualX, state ) );
		}
		state[actualNode] = false;
		actualY++;
		actualX--;
		if( ( nodeCount - actualY >= conditionNumber ) && ( actualX <= conditionNumber ) && ( actualPriceY < minimalPrice ) ) {
			nextStates.push_back( BFSState( actualNode + 1, actualPriceY, actualX, state ) );
		}
	}
	return nextStates;	
}

/*======================================Recursive function======================================*/
/*
 * Recursive function for sequential solving that finds minimal edge cut
 * actualNode - index of actual vertex
 * actualPrice - actual price of cut
 * actualX - actual count of vertices in X set
 * state - vector which describes both sets with vertices
 * solve - 	1) Add new boolean to state vector for current vertex (doesn't matter if it's true or 
 				false) and set initial values.
 			2) Compute price increments (actualPriceX for vertex in X, actualPriceY for vertex in Y).
 			3) If it's final vertex, then check if the price is less then current minimal price
 				for both option (vertex in X and Y). Two threads could modify this value, so there
 				has to be critical section handling.
 			4) If it's not final vertex try to add vertex to X and then to Y and call solve function
 				for next vertex.

 */

void solve( int actualNode, double actualPrice, int actualX, vector<bool> & state ) {
	//Initialization
	state.push_back( false );
	int actualY = actualNode - actualX;
	double actualPriceX = actualPrice;
	double actualPriceY = actualPrice;
	//Price computing
	for ( int i = 0 ; i < ( int ) graph.graph[actualNode].size( ) ; i++ ){
		if( state[graph.graph[actualNode][i].getStartNode( )] ){
			actualPriceY += graph.graph[actualNode][i].getPrice( );
		}
		else {
			actualPriceX += graph.graph[actualNode][i].getPrice( );
		}
	}
	//Adding to sets
	if( actualNode == ( nodeCount - 1 ) ) {
		if ( actualX + 1 == conditionNumber ) {
			#pragma omp critical
			{
				if ( minimalPrice > actualPriceX ){
					minimalPrice = actualPriceX;
					state[actualNode] = true;
					minimalCutState = state;
				}
			}
		}
		else if ( actualX == conditionNumber ) {
			#pragma omp critical
			{
				if ( minimalPrice > actualPriceY ){
					minimalPrice = actualPriceY;
					state[actualNode] = false;
					minimalCutState = state;
				}
			}
		}
	}
	else {
		//Adding to X
		state[actualNode] = true;
		actualX++;
		if( ( nodeCount - actualY >= conditionNumber ) && ( actualX <= conditionNumber ) && ( actualPriceX < minimalPrice ) ) {
			solve( actualNode + 1, actualPriceX, actualX, state );
		}
		state[actualNode] = false;
		actualY++;
		actualX--;
		//Adding to Y
		if( ( nodeCount - actualY >= conditionNumber ) && ( actualX <= conditionNumber ) && ( actualPriceY < minimalPrice ) ) {
			solve( actualNode + 1, actualPriceY, actualX, state );
		}
	}
}

/*==================================Start of sequential solve====================================*/
/*
 * Void that start sequential solving
 * bs - concrete BFS state 
 */

void solve( BFSState bS ) {
	solve( bS.actualNode, bS.actualPrice, bS.actualX, bS.state );
}

/*===============================Solve call with data parallelism================================*/
/*
 * 1) Creating deque "states" and add init state
 * 2) Using WHILE cycle to generate enough states (ENOUGH_STATES is defined at start of this file)
 * 3) Using pragma omp parallel for data parallelism
 */
void solve( ) {
	deque<BFSState> states;
	vector<bool> emptyVector;
	states.push_back( BFSState( 0, 0, 0, emptyVector ) );
	while( ( int )states.size( ) < ENOUGH_STATES ) {
		deque<BFSState> nextStates = states.front( ).getNextBFSStates( );
		while ( nextStates.size( ) != 0 ) {
			states.push_back( nextStates.front( ) );
			nextStates.pop_front( );
		}
		states.pop_front( );
	}
	int countOfStates = states.size( );
	#pragma omp parallel for schedule ( dynamic )
		for ( int i = 0 ; i < countOfStates ; i++ ){
			solve( states[i] );
		}
}

/*=====================================Print X and Y set========================================*/
/*
 * Formal representation for sets
 * setToFormalSet - function which convert set with vertices to string in format $name={$node_1, $node_2, ... $node_n}
 * printSets - void which gets and prints sets X and Y 
 */

string setToFormalSet( set<int> setA, string name ) {
	string result = name + "={";
	int index = 0 ;
	for ( int number : setA ) {
		if ( index == 0 ) {
			index = 1 ;
			result += to_string( number );
		}
		else{
			result += ", " + to_string( number );
		}
	}
	result += "}";
	return result;
}

void printSets( ) {
	set<int> setX;
	set<int> setY;
	for ( int i = 0 ; i < nodeCount ; i++ ) {
		if ( minimalCutState[i] ) {
			setX.insert( i );
		}
		else{
			setY.insert( i );	
		}	
	}
	cout << setToFormalSet( setX, "X" ) << endl;
	cout << setToFormalSet( setY, "Y" ) << endl << endl;
}

/*=====================================READ INPUT DATA==========================================*/
/*
 * Reading input data
 * readInput - void for reading input
 */

void readInput( ){
	cin >> nodeCount;
	cin >> edgeCount;
	cin >> conditionNumber;
	graph = Graph( nodeCount, edgeCount );
}

/*====================================WRITE OUTPUT DATA=========================================*/
/* 
 * Result visualisation
 * printOutput - void which prints table with results (minimum price and time)
 */

void printOutput( Result myResult, Result referenceResult ){
	cout << "_________________________________________________________________________________________________________" << endl << endl ;
	printSets( );
	cout << "Results:" << endl;
	cout << "+=========================+=========================+=========================+=========================+" << endl ;
	cout << '|' << setw(25) << "Detail" << '|' << setw(25) << "Measured value" << '|' << setw(25) <<  "Reference value" << '|' << setw(25) <<  "Comparison" << '|' << endl;
	cout << "+=========================+=========================+=========================+=========================+" << endl ;
	cout << '|' << setw(25) << "Price" << '|' << setw(25) << myResult.getPrice( ) << '|' << setw(25) << referenceResult.getPrice( ) << '|' << setw(24) << round(referenceResult.getPrice( ) / myResult.getPrice( ) * 100)  << "%" << '|' << endl;
	cout << "+-------------------------+-------------------------+-------------------------+-------------------------+" << endl ;
	cout << '|' << setw(25) << "Time" << '|' << setw(24) << myResult.getTime( ) << "s" << '|' << setw(24) << referenceResult.getTime( ) << "s" << '|' << setw(24) << round(referenceResult.getTime( ) / myResult.getTime( ) * 100) << "%" << '|' << endl;
	cout << "+=========================+=========================+=========================+=========================+" << endl ;
}

/*=============================================MAIN=============================================*/
/*
 * Main function with whole process
 * 1) Reading input
 * 2) Finding minimal edge cut
 * 3) Measuring time
 * 4) Showing results
 */

int main( int argc, char** argv ) {
	readInput( );
	auto start = steady_clock::now();
	solve( );
	auto end = steady_clock::now();
	double time = (double) chrono::duration_cast<chrono::nanoseconds>( end - start ).count( );
	printOutput( Result( minimalPrice, 0, time / 1e9 ) , Result( atof( argv[1] ), atof( argv[2] ), atof( argv[3] ) ) );
	return 0;
}

/*==============================================================================================*/