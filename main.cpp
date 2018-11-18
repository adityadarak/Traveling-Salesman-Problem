#include<iostream>
#include "tsp.h"
#include "twoOpt.h"
using namespace std;


int main(int argc, char** argv){
	if(argc<2){
		cout<<"No input file Submitted"<<endl;
		return 0;
	}

	string in;
	in = argv[1];
	TSP tsp(in);
	tsp.readCities();
	cout<<"Number of cities: "<<tsp.n<<endl;

	cout<<"Fill Matrix"<<endl;
	tsp.fillMatrix();

	// tsp.directFillMatrix();

	cout<<"Finding MST in Graph"<<endl;
	tsp.findMST_old();

	cout<<"Finding Perfect Matching"<<endl;
	tsp.perfect_matching();

	int best = INT_MAX;
	int bestIndex;

	for(long i=0; i<tsp.n; i++){
		if (tsp.path_vals[i][1] < best) {
			bestIndex = tsp.path_vals[i][0];
			best = tsp.path_vals[i][1];
		}
	}

	// Store best path
	tsp.create_tour(bestIndex);
	tsp.make_shorter();

	cout<<"After Christofides Algorithm:"<<endl;
	cout << "Final length: " << tsp.pathLength << endl;

	cout<<"Enter the number of iterations of genetic algorithm:"<<endl;
	int g;
	cin>>g;

	Genetic genetic(&tsp, 10, g, 50, false, in);

	const clock_t begin_time = clock(); // gets time
	genetic.run(); // runs the genetic algorithm
	cout << "\n\nTime for to run the genetic algorithm: " << float(clock () - begin_time) /  CLOCKS_PER_SEC << " seconds."; // shows time in seconds
	


}