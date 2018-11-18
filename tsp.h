
#ifndef MWM_H_
#define MWM_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <utility> // pair
#include <time.h> // time
#include <stdlib.h> // srand, rand

#include "twoOpt.h"

using namespace std;

class TSP{
private:
	struct City{
		int x;
		int y;
	};
	//Input File
	string inFile;

	//List of odd Nodes
	vector<int> odds;

	int **cost;

	void getCount();

	void findOdds();

	int minKey(int key[], bool mstSet[]);

protected:


public:
	// Number of nodes
	int n;

	// euler circuit
	vector<int>circuit;

	// Store cities and coords read in from file
	vector<City>cities;

	// Full n x n cost matrix of distances between each city
	int **graph;

	// Current shortest path length
	int pathLength;


	// Adjacency list
	// Array of n dynamic arrays, each holding a list of nodes it's index is attached to
	vector<int> *adjlist;

	// n x 2 array to store length of TSP path starting at each node
	// col 0 => starting index   col 1 => path length from that node
	int **path_vals;


	// Constructor
	TSP(string input);

	// Destructor
	~TSP();

	// Calculate distance
	int get_distance(struct City c1, struct City c2);


	// Initialization functions
	void readCities();
	void fillMatrix();

	void directFillMatrix();

	// Find MST using Prim's algorithm
	void findMST_old();


	// Find perfect matching
	void perfect_matching();

	// Find best node to start euler at
	// Doesn't create tour, just checks
	int find_best_path(int);

	// Create tour starting at specified node
	void create_tour(int);

	// Private functions implemented by create_tour() and find_best_path()
	void euler (int pos, vector<int> &);
	//void euler(int);
	void make_hamilton(vector<int> &, int&);

	// Calls twoOpt function
	void make_shorter();


	// Debugging functions
	void printCities();
	void printAdjList();
	void printResult();
	void printEuler();
	void printPath();

	// Get node count
	int get_size() {return n;};

};

typedef std::pair<std::vector<int>, int> my_pair;

struct sort_pred
{
	bool operator()(const my_pair& firstElem, const my_pair& secondElem)
	{
		return firstElem.second < secondElem.second;
	}
};



class Genetic
{
private:
	string in;
	TSP* tsp; // the graph
	std::vector<pair<vector<int>, int> > population; // each element is a pair: vector and total cost
	pair<vector<int>, int> bestPop;
	int size_population; // size of population
	int real_size_population; // real size population
	int generations; // amount of generations
	int mutation_rate; // mutation rate
	bool show_population; // flag to show population
private:
	void initialPopulation(); // generates the initial population
public:
	Genetic(TSP* tsp, int amount_population, int generations, int mutation_rate, bool show_population = true, string in="out"); // constructor
	int isValidSolution(std::vector<int>& solution); // checks if a solution is valid
	void showPopulation(); // shows population
	void doCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2);
	void doNearestNeighborCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2);
	void doNewCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2);
	void doMCycleCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2);
	void mutation(vector<int>& child);
	void crossOver(std::vector<int>& parent1, std::vector<int>& parent2, vector< pair<vector<int> ,int> >& newPopulation); // makes the crossover
	void insertBinarySearch(std::vector<int>& child, int total_cost, vector< pair<vector<int> ,int> >& newPopulation); // uses binary search to insert
	void run(); // runs genetic algorithm
	int getCostBestSolution(); // returns cost of the best solution
	bool existsChromosome(const std::vector<int> & v); // checks if exists the chromosome
};



#endif