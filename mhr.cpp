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
#include <omp.h>
#include <mpi.h>
#include <utility>
#include <bits/stdc++.h> 

using namespace std;
using namespace chrono;

/*==========================================CONSTANTS===========================================*/
/*
 * MAX_THREADS - count of threads 
 * ENOUGH_STATES_MASTER - count of states, that master must have in container
 * ENOUGH_STATES_SLAVE - count of states, that slave must have in container to start parallel run
 */

#define MAX_THREADS 8
#define ENOUGH_STATES_MASTER 80
#define ENOUGH_STATES_SLAVE 800
/*================================GLOBAL VARIABLE - READ & WRITE================================*/
/*
 * minimalPrice - global variable to store the price of minimal edge cut
 * minimalCutState - global vector which describes both sets with vertices
 */

double minimalPrice = 0.0;
vector<bool> minimalCutState;

/*================================Convert state for MPI functions===============================*/
/*
 * Funtion to convert vector to double for MPI sending
 */

double vectorToDouble( vector<bool> state ) {
    double result = 1;
    for ( int i = 0 ; i < (int)state.size( ) ; i++ ) {
        if ( state[i] ) {
            result = result * 2 + 1;
        }
        else {
            result *= 2;
        }
    }
    return result; 
}

/*----------------------------------------------------------------------------------------------*/
/*
 * Funtion to convert double to vector from MPI sending
 */

vector<bool> doubleToVector( double stateNum ) {
    vector<bool> result;
    while ( stateNum > 1.5 ) {
        if ( fmod( stateNum, 2 ) > 0.5 ) {
            result.push_back( true );
            stateNum = ( stateNum - 1 ) / 2 ;
        }
        else {
            result.push_back( false );
            stateNum = stateNum / 2 ;
        }
    }
    reverse( result.begin( ), result.end( ) );
    return result;
}

/*=====================================Some input functions=====================================*/
/*
 * Function for string replacement
 */

string ReplaceAll(string str, const string& from, const string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
    return str;
}

/*----------------------------------------------------------------------------------------------*/
/*
 * Function that returns parts separated by a space
 * Parts of string line are in deque
 */

deque<string> splitLineFromString( string line ) {
	deque<string> result ;
	bool space_found = true ;
	while( space_found ) {
		space_found = false ;
		for( int i = 0 ; i <= (int)line.size( ) ; i++ ) {
			if ( line[i] == ' ' ) {
				space_found = true ;
				result.push_back(line.substr(0,i));
				line = line.substr(i+1);
				break ;
			}
		}
	}
	result.push_back( line );
    return result;
}  

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
	Graph( int nodes, int edges, char* filename );
	void addEdge( Edge edge );
	vector<Edge> getEdgesToNode( int node );
	void loadGraph( char* filename );
	vector<vector<Edge>> graph;
private:
	int nodeCount;
	int edgeCount;
};

/*----------------------------------------------------------------------------------------------*/

Graph::Graph( int nodes, int edges, char* filename ) {
	for ( int i = 0 ; i < nodes ; i++ ){
		vector<Edge> vectorTmp;
		this->graph.push_back( vectorTmp );
	}
	this->nodeCount = nodes;
	this->edgeCount = edges;
	this->loadGraph( filename );
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

void Graph::loadGraph( char* filename ) {
	if ( ( this->nodeCount != 0 ) && ( this->edgeCount != 0 ) ) {
		ifstream infile( filename );
		string line;
		getline( infile, line );
		int startNode;
		int endNode;
		double price;
		for ( int i = 0 ; i < nodeCount * edgeCount / 2 ; i++ ) {
        	getline( infile, line );
        	line = ReplaceAll(line, "    ", " ");
        	line = ReplaceAll(line, "   ", " ");
        	line = ReplaceAll(line, "  ", " ");
        	deque<string> values = splitLineFromString( line );
			startNode = atoi( values[1].c_str( ) );
			endNode = atoi( values[2].c_str( ) );
			price = atof( values[3].c_str( ) );
			minimalPrice += price;
			this->graph[endNode].push_back( Edge( startNode, endNode, price ) );
		}
		minimalPrice += 1.0;
	}
}

/*=================================GLOBAL VARIABLE - READ ONLY==================================*/
/* 
 * graph - instance of Graph class which has information about graph (vertices and edges)
 * nodeCount - count of vertices
 * edgeCount - count of edges
 * conditionNumber - count of vertices that the set X must contain
 */

char *emptyCharPointer;
Graph graph = Graph( 0, 0, emptyCharPointer );
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
	pair<deque<BFSState>,bool> getNextBFSStates( );
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

pair<deque<BFSState>,bool> BFSState::getNextBFSStates( ) {
	deque<BFSState> nextStates;
	bool flag = false;
	if ( this->actualNode == ( nodeCount - 1 ) ) {
		nextStates.push_back( *this );
		flag = true;
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
	return make_pair(nextStates, flag);	
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
	if ( (int)state.size( ) < nodeCount ) {
		state.push_back( false );
	}
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


/*==========================================MPI part=============================================*/
/*
 * Struct that represents bfs state
 * For sending state from master to slave in MPI 
 */

typedef struct bfs_s{
	int actualNode; 
	double actualPrice; 
	int actualX; 
	double state;
}bfs_state;

/*----------------------------------------------------------------------------------------------*/
/*
 * Struct that represents result of solve function
 * For sending result from slave back to master in MPI 
 */

typedef struct s_r{
	double price;
	double state;
}solve_result;

/*----------------------------------------------------------------------------------------------*/
/*
 * MPI solve function
 * Comment in code
 */

void solveMPI( int argc, char** argv  ) {
	//Get rank of process and count of processes
	int my_rank, num_procs;
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size( MPI_COMM_WORLD, &num_procs);
	//Creating tags for communication
	static const int tag_done = 1;
	static const int tag_finished = 2;
	static const int tag_more_work = 3;
	static const int tag_work = 4;
	//Creating MPI datatypes, which are sending in communication
	const int BFSItemsCount = 4;
	const int ResultItemsCount = 2;
	int BFSBlocklengths[] = { 1, 1, 1, 1};
	int ResultBlocklengths[] = { 1, 1 };
	MPI_Aint BFSOffsets[4];
	MPI_Aint ResultOffsets[2];
	BFSOffsets[0] = offsetof( bfs_state, actualNode );
	BFSOffsets[1] = offsetof( bfs_state, actualPrice );
	BFSOffsets[2] = offsetof( bfs_state, actualX );
	BFSOffsets[3] = offsetof( bfs_state, state );
	ResultOffsets[0] = offsetof( solve_result, price );
	ResultOffsets[1] = offsetof( solve_result, state );
	MPI_Datatype BFStypes[] = { MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE };
	MPI_Datatype Resulttypes[] = { MPI_DOUBLE, MPI_DOUBLE };
	MPI_Datatype mpi_bfs_type;
	MPI_Datatype mpi_result_type;
	MPI_Type_create_struct( BFSItemsCount, BFSBlocklengths, BFSOffsets, BFStypes, &mpi_bfs_type );
	MPI_Type_commit( &mpi_bfs_type );
	MPI_Type_create_struct( ResultItemsCount, ResultBlocklengths, ResultOffsets, Resulttypes, &mpi_result_type );
	MPI_Type_commit( &mpi_result_type );
	// If rank of process is zero, then the process is MASTER
	// If rank of process is not zero, then the process is SLAVE
	if ( my_rank == 0 ) {
		// MASTER process
		auto start = steady_clock::now();
		// First fill the container with states
		deque<BFSState> states;
		vector<bool> emptyVector;
		states.push_back( BFSState( 0, 0, 0, emptyVector ) );
		int flagCounter = 0;
		// If deque has enough states, or if there are no updates, or if the deque is empty, then end filling
		// Else give next states to deque
		while( ( ( int )states.size( ) < ENOUGH_STATES_MASTER ) && ( flagCounter < ENOUGH_STATES_SLAVE ) && ( states.size( ) != 0 ) ) {
			pair<deque<BFSState>,bool> deboPair = states.front( ).getNextBFSStates( );
			deque<BFSState> nextStates = deboPair.first;
			if ( deboPair.second ) {
				flagCounter++;
			}
			else {
				flagCounter = 0;
			}
			while ( nextStates.size( ) != 0 ) {
				states.push_back( nextStates.front( ) );
				nextStates.pop_front( );
			}
			states.pop_front( );
		}
		// Send info, that there is work to all SLAVE processes
		// Send work to all SLAVEs 
		for( int dest = 1; dest < num_procs; dest++ ) {
			int msg = 1;
			MPI_Send( &msg, 1, MPI_INT, dest, tag_more_work, MPI_COMM_WORLD );
			BFSState bfsState = states.front( );
			states.pop_front( );
			bfs_state sendItem;
			sendItem.actualNode = bfsState.actualNode;
			sendItem.actualPrice = bfsState.actualPrice;
			sendItem.actualX = bfsState.actualX;
			sendItem.state = vectorToDouble( bfsState.state );
			MPI_Send( &sendItem, 1, mpi_bfs_type, dest, tag_work, MPI_COMM_WORLD );
		}
		// Variable for store count of active processes
		int workingSlaves = num_procs - 1;
		// If there is some working SLAVE process, then send him work, or info that there is no work
		while( workingSlaves != 0 ) {
			MPI_Status status;
			solve_result result;
			// Recieve result from SLAVE that has completed its work
			MPI_Recv( &result, 1, mpi_result_type, MPI_ANY_SOURCE, tag_done, MPI_COMM_WORLD, &status );
			// If the result is better then actual minimal price, then update this price
			if( result.price < minimalPrice ){
				minimalPrice = result.price;
				minimalCutState = doubleToVector( result.state );
			}
			// If there is no work, send tag_finished to SLAVE and update count of active processes
			if( states.size( ) == 0 ) {
				int msg = 1;
				MPI_Send( &msg, 1, MPI_INT, status.MPI_SOURCE, tag_finished, MPI_COMM_WORLD );
				workingSlaves--;
			}
			// Else send to SLAVE next work
			else {
				int msg = 1;
				MPI_Send( &msg, 1, MPI_INT, status.MPI_SOURCE, tag_more_work, MPI_COMM_WORLD );
				BFSState bfsState = states.front( );
				states.pop_front( );
				bfs_state sendItem;
				sendItem.actualNode = bfsState.actualNode;
				sendItem.actualPrice = bfsState.actualPrice;
				sendItem.actualX = bfsState.actualX;
				sendItem.state = vectorToDouble( bfsState.state );
				MPI_Send( &sendItem, 1, mpi_bfs_type, status.MPI_SOURCE, tag_work, MPI_COMM_WORLD );
			}
		}
		// At the end print the result of solving
		auto end = steady_clock::now();
		double time = (double) chrono::duration_cast<chrono::nanoseconds>( end - start ).count( );
		printOutput( Result( minimalPrice, 0, time / 1e9 ) , Result( atof( argv[1] ), atof( argv[2] ), atof( argv[3] ) ) );
	}
	else {
		// SLAVE process
		while( true ){
			MPI_Status status;
			int msg;
			// First recieve info if there is work
			// If there is no work, then end
			MPI_Recv( &msg, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			if( status.MPI_TAG == tag_finished ){
				return;
			}
			// Recieve work from MASTER
			bfs_state bfsStateElement;
			MPI_Recv( &bfsStateElement, 1, mpi_bfs_type, 0, tag_work, MPI_COMM_WORLD, &status );
			// Fill deque with states for parallel run
			deque<BFSState> states;
			vector<bool> emptyVector;
			states.push_back( BFSState(bfsStateElement.actualNode, bfsStateElement.actualPrice, bfsStateElement.actualX, doubleToVector( bfsStateElement.state ) ) );
			int flagCounter = 0;
			// If deque has enough states, or if there are no updates, or if the deque is empty, then end filling
			// Else give next states to deque
			while( (  ( int )states.size( ) < ENOUGH_STATES_SLAVE ) && ( flagCounter < ENOUGH_STATES_SLAVE ) && ( states.size( ) != 0 ) ) {
				pair<deque<BFSState>,bool> deboPair = states.front( ).getNextBFSStates( );
				deque<BFSState> nextStates = deboPair.first;
				if ( deboPair.second ) {
					flagCounter++;
				}
				else {
					flagCounter = 0;
				}
				while ( nextStates.size( ) != 0 ) {
					states.push_back( nextStates.front( ) );
					nextStates.pop_front( );
				}
				states.pop_front( );
			}
			// Parallel solving			
			int countOfStates = states.size( );
			#pragma omp parallel for schedule ( dynamic ) num_threads( MAX_THREADS )
				for ( int i = 0 ; i < countOfStates ; i++ ){
					solve( states[i] );	
				}

			// Send result to MASTER process
			solve_result result;
			result.price = minimalPrice;
			result.state = vectorToDouble( minimalCutState );
			MPI_Send( &result, 1, mpi_result_type, 0, tag_done, MPI_COMM_WORLD );
		}
	}
}

/*=====================================READ INPUT DATA==========================================*/
/*
 * Reading input data
 * readInput - void for reading input
 */

void readInput( char* filename ){
    ifstream infile( filename );
	string line;
	getline( infile, line );
    deque<string> values = splitLineFromString( line );
    nodeCount = atoi( values[0].c_str( ) );
	edgeCount = atoi( values[1].c_str( ) ) ;
	conditionNumber = atoi( values[2].c_str () );
	graph = Graph( nodeCount, edgeCount, filename );
}

/*=============================================MAIN=============================================*/
/*
 * MPI initialization
 * All processes read the input and make their own graph
 * Then all processes start solveMPI void
 * MPI finalization 
 */

int main( int argc, char** argv ) {
	MPI_Init( &argc, &argv );
	readInput( argv[4] );
	solveMPI( argc, argv );
	MPI_Finalize( );
	return 0;
}

/*==============================================================================================*/