#include <iostream>
#include <string.h>
#include <cmath>
#include <unordered_map>
#include <set>
#include <string>
#include <iomanip>
#include <vector>
#include <chrono>

using namespace std;
using namespace chrono;

/*================================GLOBAL VARIABLE - READ & WRITE================================*/

double minimalPrice = 0.0;
vector<bool> minimalCutState;
vector<bool> state;
double calling = 0.0 ;

/*=========================================RESULT class=========================================*/

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

class Graph{
public:
	Graph( int nodes, int edges );
	void addEdge( Edge edge );
	vector<Edge> getEdgesToNode( int node );
	vector<vector<Edge>> getGraph( );
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

vector<vector<Edge>> Graph::getGraph( ) {
	return this->graph;
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

Graph graph = Graph( 0, 0 );
int nodeCount = 0;
int edgeCount = 0;
int conditionNumber = 0;

/*======================================Recursion function======================================*/

void solve( int actualNode, double actualPrice, int actualX ) {
	//Počáteční inicializace
	calling++ ;
	int actualY = actualNode - actualX;
	double actualPriceX = actualPrice;
	double actualPriceY = actualPrice;
	//Napočítání cen
	for ( int i = 0 ; i < ( int ) graph.graph[actualNode].size( ) ; i++ ){
		if( state[graph.graph[actualNode][i].getStartNode( )] ){
			actualPriceY += graph.graph[actualNode][i].getPrice( );
		}
		else {
			actualPriceX += graph.graph[actualNode][i].getPrice( );
		}
	}
	//Přidání
	if( actualNode == ( nodeCount - 1 ) ) {
		if ( actualX + 1 == conditionNumber ) {
			if ( minimalPrice > actualPriceX ){
				minimalPrice = actualPriceX;
				state[actualNode] = true;
				minimalCutState = state;
			}
		}
		else if ( actualX == conditionNumber ) {
			if ( minimalPrice > actualPriceY ){
				minimalPrice = actualPriceY;
				state[actualNode] = false;
				minimalCutState = state;
			}
		}
	}
	else {
		//Přidávám do X
		state[actualNode] = true;
		actualX++;
		if( ( nodeCount - actualY >= conditionNumber ) && ( actualX <= conditionNumber ) && ( actualPriceX < minimalPrice ) ) {
			solve( actualNode + 1, actualPriceX, actualX );
		}
		state[actualNode] = false;
		actualY++;
		actualX--;
		//Přidávám do Y
		if( ( nodeCount - actualY >= conditionNumber ) && ( actualX <= conditionNumber ) && ( actualPriceY < minimalPrice ) ) {
			solve( actualNode + 1, actualPriceY, actualX );
		}
	}
}

/*=====================================Initial solve call=======================================*/

void solve( ) {
	vector<bool> falseVector;
	for ( int i = 1 ; i <= nodeCount ; i++ ) {
		falseVector.push_back( false );
	}
	state = falseVector;
	solve( 0, 0, 0 );
}

/*=====================================Print X and Y set========================================*/

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

void readInput( ){
	cin >> nodeCount;
	cin >> edgeCount;
	cin >> conditionNumber;
	graph = Graph( nodeCount, edgeCount );
}

/*====================================WRITE OUTPUT DATA=========================================*/

void printOutput( Result myResult, Result referenceResult ){
	cout << "_________________________________________________________________________________________________________" << endl << endl ;
	printSets( );
	cout << "Results:" << endl;
	cout << "+=========================+=========================+=========================+=========================+" << endl ;
	cout << '|' << setw(25) << "Detail" << '|' << setw(25) << "Measured value" << '|' << setw(25) <<  "Reference value" << '|' << setw(25) <<  "Comparison" << '|' << endl;
	cout << "+=========================+=========================+=========================+=========================+" << endl ;
	cout << '|' << setw(25) << "Price" << '|' << setw(25) << myResult.getPrice( ) << '|' << setw(25) << referenceResult.getPrice( ) << '|' << setw(24) << round(referenceResult.getPrice( ) / myResult.getPrice( ) * 100)  << "%" << '|' << endl;
	cout << "+-------------------------+-------------------------+-------------------------+-------------------------+" << endl ;
	cout << '|' << setw(25) << "Complexity" << '|' << setw(25) << myResult.getComplexity( ) << '|' << setw(25) << referenceResult.getComplexity( ) << '|' << setw(24) << round(referenceResult.getComplexity( ) / myResult.getComplexity( ) * 100)  << "%" << '|' << endl;
	cout << "+-------------------------+-------------------------+-------------------------+-------------------------+" << endl ;
	cout << '|' << setw(25) << "Time" << '|' << setw(24) << myResult.getTime( ) << "s" << '|' << setw(24) << referenceResult.getTime( ) << "s" << '|' << setw(24) << round(referenceResult.getTime( ) / myResult.getTime( ) * 100) << "%" << '|' << endl;
	cout << "+=========================+=========================+=========================+=========================+" << endl ;
}

/*=============================================MAIN=============================================*/

int main( int argc, char** argv ) {
	auto start = steady_clock::now();
	readInput( );
	solve( );
	auto end = steady_clock::now();
	double time = (double) chrono::duration_cast<chrono::nanoseconds>( end - start ).count( );
	printOutput( Result( minimalPrice, calling, time / 1000000000 ) , Result( atof( argv[1] ), atof( argv[2] ), atof( argv[3] ) ) );
	return 0;
}

/*==============================================================================================*/