
#include "tsp.h"

TSP::TSP(string input){
	inFile = input;

	getCount();

	graph = new int*[n];
	for (int i = 0; i < n; i++) {
		graph[i] = new int[n];
		for (int j = 0; j < n; j++) graph[i][j] = 0;
	}

	cost = new int*[n];
	for (int i = 0; i < n; i++) {
		cost[i] = new int[n];
	}

	path_vals = new int*[n];
	for (int i = 0; i < n; i++) {
		path_vals[i] = new int[n];
	}

	// Adjacency list
	adjlist = new vector<int> [n];
}

TSP::~TSP(){
	for (int i = 0; i < n; i++) {
		delete [] graph[i];
		delete [] cost[i];
		delete [] path_vals[i];
	}
	delete [] path_vals;
	delete [] graph;
	delete [] cost;
	delete [] adjlist;
}

void TSP::getCount(){
	int count = 0;
	cout<<"Getting Count"<<endl;
	ifstream inStream(inFile.c_str());
	// inStream.open(inFile.c_str() ios::in);

	if (!inStream) {
	  cerr << "Can't open input file " << inFile << endl;
	  exit(1);
	}
	std::string unused;
	while ( std::getline(inStream, unused) )
	   ++count;
	n = count;
	cout<<"No. of cities: "<<n<<endl;
	inStream.close();
};

void TSP::readCities(){
	ifstream inStream(inFile.c_str());
	cout<<"Reading Cities"<<endl;
	// inStream.open(inFile.c_str(), ios::in);
	if (!inStream) {
	  cerr << "Can't open input file " << inFile << endl;
	  exit(1);
	}
	int c, x, y;
	int i = 0;
	while (!inStream.eof() ) {
		inStream >> c >> x >> y;
		// cout<<c<<" "<<x<<" "<<y<<endl;
		// Push back new city to vector
		struct City c = {x, y};
		cities.push_back(c);
		i++;
	}
	inStream.close();
};

int TSP::get_distance(struct TSP::City c1, struct TSP::City c2) {
	int dx = pow((float)(c1.x - c2.x), 2);
	int dy = pow((float)(c1.y - c2.y), 2);
	return (floor((float) (sqrt(dx + dy)) + 0.5));
};

void TSP::fillMatrix(){
	for(int i=0;i<n;i++){
		for(int j=i; j<n;j++){
			graph[i][j] = graph[j][i] = get_distance(cities[i], cities[j]);
		}
	}
};

void TSP::directFillMatrix(){
	ifstream inStream(inFile.c_str());
	cout<<"Reading Cities"<<endl;
	// inStream.open(inFile.c_str(), ios::in);
	if (!inStream) {
	  cerr << "Can't open input file " << inFile << endl;
	  exit(1);
	}

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			inStream>>graph[i][j];
		}
	}
}

void TSP::findMST_old() {
	int key[n];   // Key values used to pick minimum weight edge in cut
	bool in_mst[n];  // To represent set of vertices not yet included in MST
	int parent[n];

	// For each vertex v in V
	for (int v = 0; v < n; v++) {
		// Initialize all keys to infinity
		key[v] = std::numeric_limits<int>::max();

		// Mark as not being in mst yet
		in_mst[v] = false;
	}

	// Node 0 is the root node so give it the lowest distance (key)
	key[0] = 0;
	parent[0] = -1; // First node is always root of MST

	for (int i = 0; i < n - 1; i++) {
		// Find closest remaining (not in tree) vertex
		// TO DO : This would be better represented by heap/pqueue
		int v = minKey(key, in_mst);

		// Add vertex v to the MST
		in_mst[v] = true;

		// Look at each vertex u adjacent to v that's not yet in mst
		for (int u = 0; u < n; u++) {
			if (graph[v][u] && in_mst[u] == false && graph[v][u] < key[u]) {
				// Update parent index of u
				parent[u] = v;

				// Update the key only if dist is smaller than key[u]
				key[u] = graph[v][u];
			}
		}
	}

	// map relations from parent array onto matrix
	for (int v1 = 0; v1 < n; v1++) {
		// there is an edge between v1 and parent[v1]
		int v2 = parent[v1];
		if (v2 != -1) {
			adjlist[v1].push_back(v2);
			adjlist[v2].push_back(v1);
		}
	}
};

// findMST helper function
int TSP::minKey(int key[], bool mstSet[]) {
	// Initialize min value
	int min = std::numeric_limits<int>::max();
	int min_index;
	for (int v = 0; v < n; v++)
		if (mstSet[v] == false && key[v] < min) {
			min = key[v];
			min_index = v;
		}
	return min_index;
};

void TSP::findOdds() {
	/////////////////////////////////////////////////////
	// Find nodes with odd degrees in T to get subgraph O
	/////////////////////////////////////////////////////

	// store odds in new vector for now
	for (int r = 0; r < n; r++) {
		//cities[r].isOdd = ((adjlist[r].size() % 2) == 0) ? 0 : 1;
		if ((adjlist[r].size() % 2) != 0 ) {
			odds.push_back(r);
		}
	}
}

void TSP::perfect_matching() {
	int closest, length; //int d;
	std::vector<int>::iterator tmp, first;

	// Find nodes with odd degrees in T to get subgraph O
	findOdds();

	// for each odd node
	while (!odds.empty()) {
		first = odds.begin();
		vector<int>::iterator it = odds.begin() + 1;
		vector<int>::iterator end = odds.end();
		length = std::numeric_limits<int>::max();
		for (; it != end; ++it) {
			// if this node is closer than the current closest, update closest and length
			if (graph[*first][*it] < length) {
				length = graph[*first][*it];
				closest = *it;
				tmp = it;
			}
		}	// two nodes are matched, end of list reached
		adjlist[*first].push_back(closest);
		adjlist[closest].push_back(*first);
		odds.erase(tmp);
		odds.erase(first);
	}
}


// Take reference to a path vector
// so can either modify actual euler path or a copy of it
void TSP::euler (int pos, vector<int> &path) {
	/////////////////////////////////////////////////////////
	// Based on this algorithm:
	//	http://www.graph-magics.com/articles/euler.php
	// we know graph has 0 odd vertices, so start at any vertex
	// O(V+E) complexity
	/////////////////////////////////////////////////////////

	// make copy of original adjlist to use/modify
	vector<int> *temp = new vector<int> [n];
	for (int i = 0; i < n; i++) {
		temp[i].resize(adjlist[i].size());
		temp[i] = adjlist[i];
	}

	path.clear();

	// Repeat until the current vertex has no more neighbors and the stack is empty.
	stack<int> stk;
	while (!stk.empty() || temp[pos].size() > 0 ) {
		// If current vertex has no neighbors -
		if (temp[pos].size() == 0) {
			// add it to circuit,
			path.push_back(pos);
			// remove the last vertex from the stack and set it as the current one.
			int last = stk.top();
			stk.pop();
			pos = last;
		}
		// Otherwise (in case it has neighbors)
		else {
			// add the vertex to the stack,
			stk.push(pos);
			// take any of its neighbors,
			int neighbor = temp[pos].back();
			// remove the edge between selected neighbor and that vertex,
			temp[pos].pop_back();
	        for (unsigned int i = 0; i < temp[neighbor].size(); i++)
	            if (temp[neighbor][i] == pos) { // find position of neighbor in list
	        	    temp[neighbor].erase (temp[neighbor].begin() + i); // remove it
	                break;
	            }
			// and set that neighbor as the current vertex.
	        pos = neighbor;
		}
	}
	path.push_back(pos);
}


void TSP::make_hamilton(vector<int> &path, int &path_dist) {
	// remove visited nodes from Euler tour
	bool visited[n]; // boolean value for each node if it has been visited yet
	memset(visited, 0, n * sizeof(bool));

	path_dist = 0;

	int root = path.front();
	vector<int>::iterator curr = path.begin();
	vector<int>::iterator next = path.begin()+1;
	visited[root] = true;

	// loop until the end of the circuit list is reached
	while ( next != path.end() ) {
		// if we haven't been to the next city yet, go there
		if (!visited[*next]) {
			path_dist += graph[*curr][*next];
			curr = next;
			visited[*curr] = true;
			next = curr + 1;
		}else {
			next = path.erase(next); // remove it
		}
	}
	// cout<<path_dist<<endl;
	// add the distance back to the root
	path_dist += graph[*curr][*next];
}





void TSP::create_tour(int pos){
	// call euler with actual circuit vector
	euler(pos, circuit);

	// make it hamiltonian
	// pass actual vars
	make_hamilton(circuit, pathLength);
	// cout<<"Pathlength after Christofides: "<<pathLength<<endl;
}


// Does euler and hamilton but doesn't modify any variables
// Just finds path length from the node specified and returns it
int TSP::find_best_path (int pos) {

	// create new vector to pass to euler function
	vector<int>path;
	euler(pos, path);

	// make it hamiltonian, pass copy of vars
	int length;
	make_hamilton(path, length);

	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);



	return length;
}

void TSP::make_shorter(){
	// Modify circuit & pathLength
	twoOpt(graph, circuit, pathLength, n);
}


Genetic::Genetic(TSP* tsp, int size_population, int generations, int mutation_rate, bool show_population, string in)
{
	if(size_population < 1) // checks if size of population is less than 1
	{
		cout << "Error: size_population < 1\n";
		exit(1);
	}
	else if(mutation_rate < 0 || mutation_rate > 100) // checks if mutation rate is less than 0
	{
		cout << "Error: mutation_rate must be >= 0 and <= 100\n";
		exit(1);
	}
		
	this->in = in;
	this->tsp = tsp;
	this->size_population = size_population;
	this->real_size_population = 0;
	this->generations = generations;
	this->mutation_rate = mutation_rate;
	this->show_population = show_population;
}



// checks if is a valid solution, then return total cost of path else return -1
int Genetic::isValidSolution(vector<int>& solution)
{
	int total_cost = 0;
	set<int> set_solution;
	
	// checks if not contains elements repeated
	for(int i = 0; i < tsp->n; i++)
		set_solution.insert(solution[i]);
	
	if(set_solution.size() != (unsigned)tsp->n)
		return -1;

	// checks if connections are valid
	for(int i = 0; i < tsp->n; i++)
	{
		if(i + 1 <  tsp->n)
		{
			if(solution[i]>=tsp->n || solution[i]<0){
				cout<<solution[i]<<endl;
			}
			if(solution[i+1]>=tsp->n  || solution[i+1]<0){
				cout<<solution[i]<<endl;
			}
			int cost = tsp->graph[solution[i]][solution[i+1]];

			// checks if exists connection
			if(cost == -1)
				return -1;
			else
				total_cost += cost;
		}
		else
		{
			// int cost = tsp->graph[solution[i]][solution[0]];
			
			// // checks if exists connection
			// if(cost == -1)
			// 	return -1;
			// else
			// 	total_cost += cost;
			break;
		}
	}
	return total_cost;
}


bool Genetic::existsChromosome(const vector<int> & v)
{
	// checks if exists in the population
	for(vector<pair<vector<int>, int> >::iterator it=population.begin(); it!=population.end(); ++it)
	{
		const vector<int>& vec = (*it).first; // gets the vector
		if(equal(v.begin(), v.end(), vec.begin())) // compares vectors
			return true;
	}
	return false;
}


void Genetic::initialPopulation() // generates the initial population
{
	vector<int> parent;
	
	// inserts initial vertex in the parent
	parent = tsp->circuit;


		
	int total_cost = isValidSolution(parent);

	bestPop = make_pair(parent, total_cost);
	
	if(total_cost != -1) // checks if the parent is valid
	{
		population.push_back(make_pair(parent, total_cost)); // inserts in the population
		real_size_population++; // increments real_size_population
	}
	
	// makes random permutations "generations" times
	for(int i = 0; i < generations; i++)
	{
		// generates a random permutation
		random_shuffle(parent.begin() + 1, parent.begin() + (rand() % (tsp->n - 1) + 1));
		
		int total_cost = isValidSolution(parent); // checks if solution is valid
		
		// checks if permutation is a valid solution and if not exists
		if(total_cost != -1 && !existsChromosome(parent))
		{
			population.push_back(make_pair(parent, total_cost)); // add in population
			real_size_population++; // increments real_size_population in the unit
		}
		if(real_size_population == size_population) // checks size population
			break;
	}
	
	// checks if real_size_population is 0
	if(real_size_population == 0)
		cout << "\nEmpty initial population";
	else
		sort(population.begin(), population.end(), sort_pred()); // sort population
	cout<<"Initial Population Generated"<<endl;
}


void Genetic::showPopulation()
{
	cout << "\nShowing solutions...\n\n";
	for(vector<pair<vector<int>, int> >::iterator it=population.begin(); it!=population.end(); ++it)
	{
		const vector<int>& vec = (*it).first; // gets the vector
		
		for(int i = 0; i < tsp->n; i++)
			cout << vec[i] << " ";
		cout << tsp->circuit[0];
		cout << " | Cost: " << (*it).second << "\n\n";
	}
	cout << "\nPopulation size: " << real_size_population << endl;
}


// inserts in the vector using binary search
void Genetic::insertBinarySearch(vector<int>& child, int total_cost, vector< pair<vector<int> ,int> >& newPopulation)
{
	int imin = 0;
	int imax = newPopulation.size()-1;
	
	while(imax >= imin)
	{
		int imid = imin + (imax - imin) / 2;
		
		if(total_cost == newPopulation[imid].second)
		{
			newPopulation.insert(newPopulation.begin() + imid, make_pair(child, total_cost));
			return;
		}
		else if(total_cost > newPopulation[imid].second)
			imin = imid + 1;
		else
			imax = imid - 1;
	}
	newPopulation.insert(newPopulation.begin() + imin, make_pair(child, total_cost));
}


/*
	Makes the crossover
	This crossover selects two random points
	These points generates substrings in both parents
	The substring inverted of parent1 is placed in parent2 and
	the substring inverted of parent2 is placed in parent1
	
	Example:
		parent1: 1 2 3 4 5
		parent2: 1 2 4 5 3
		
		substring in parent1: 2 3 4
		substring in parent2: 2 4 5
		
		substring inverted in parent1: 4 3 2
		substring inverted in parent2: 5 4 2
		
		child1: 1 5 4 2 5
		child2: 1 4 3 2 3
		
		Children are invalids: 5 appears 2x in child1 and 3 appears 2x in child2
		Solution: map of genes that checks if genes are not used
*/
void Genetic::doCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2){
	map<int, int> genes1, genes2;
	
	for(int i = 0; i < tsp->n; i++)
	{
		// initially the genes not are used
		genes1[parent1[i]] = 0;
		genes2[parent2[i]] = 0;
	}
	
	// generates random points
	int point1 = rand() % (tsp->n - 1) + 1;
	int point2 = rand() % (tsp->n - point1) + point1;
	
	// adjusts the points if they are equal
	if(point1 == point2)
	{
		if(point1 - 1 > 1)
			point1--;
		else if(point2 + 1 < tsp->n)
			point2++;
		else
		{
			// point1 or point2 ?? random...
			int point = rand() % 10 + 1; // number in the range 1 to 10
			if(point <= 5)
				point1--;
			else
				point2++;
		}
	}
	
	// generates childs
	
	// until point1, child1 receives genes of the parent1
	// and child2 receives genes of the parent2
	for(int i = 0; i < point1; i++)
	{
		// adds genes
		child1.push_back(parent1[i]);
		child2.push_back(parent2[i]);
		// marks genes
		genes1[parent1[i]] = 1;
		genes2[parent2[i]] = 1;
	}
	
	// marks remaining genes
	for(int i = point2 + 1; i < tsp->n; i++)
	{
		genes1[parent1[i]] = 1;
		genes2[parent2[i]] = 1;
	}
		
	// here is the substring inverted
	// child1 receives genes of the parent2 and
	// child2 receives genes of the parent1
	for(int i = point2; i >= point1; i--)
	{
		if(genes1[parent2[i]] == 0) // if the gene is not used
		{
			child1.push_back(parent2[i]);
			genes1[parent2[i]] = 1; // marks the gene	
		}
		else
		{
			// if the gene already is used, chooses gene that is not used
			for(map<int, int>::iterator it = genes1.begin(); it != genes1.end(); ++it)
			{
				if(it->second == 0) // checks if is not used
				{
					child1.push_back(it->first);
					genes1[it->first] = 1; // marks as used
					break; // left the loop
				}
			}
		}
		
		if(genes2[parent1[i]] == 0) // if the gene is not used
		{
			child2.push_back(parent1[i]);
			genes2[parent1[i]] = 1; // marks the gene
		}
		else
		{
			// if the gene already is used, chooses gene that is not used
			for(map<int, int>::iterator it = genes2.begin(); it != genes2.end(); ++it)
			{
				if(it->second == 0) // checks if is not used
				{
					child2.push_back(it->first);
					genes2[it->first] = 1; // marks as used
					break; // left the loop
				}
			}
		}
	}
	
	// ramaining genes: child1 receives genes of the parent1
	// and child2 receives genes of the parent2
	for(int i = point2 + 1; i < tsp->n; i++)
	{
		child1.push_back(parent1[i]);
		child2.push_back(parent2[i]);
	}
}

void Genetic::doNearestNeighborCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2){
	int c = rand()%tsp->n;
	// cout<<"c:"<<c<<endl;
	int q =c;
	vector<int> p1 =parent1, p2=parent2;
	child1.push_back(c);
	while(p1.size()!=1){
		vector<int>::iterator it1;
		it1 = p1.begin();
		while(*it1 != c){
			it1++;
		}
		vector<int>::iterator it2;
		it2 = p2.begin();
		while(*it2 != c){
			it2++;
		}
		// cout<<it1<<" "<<it2<<endl;
		
		int dx, dy;
		if(it1+1!=p1.end()){
			dx = *(it1+1);
		}
		else{
			dx = p1[0]; 
		}
		if(it2+1!=p2.end()){
			dy = *(it2+1);
		}
		else{
			dy = p2[0]; 
		}

		// int dx = p1[(*(it1+1))%p1.size()];
		// int dy = p2[(*(it2+1))%p1.size()];
		// cout<<dx<<" "<<dy<<endl;
		
        c = tsp->graph[c][dx] > tsp->graph[c][dy] ? dy:dx;
        p1.erase(it1);
		p2.erase(it2);
        child1.push_back(c);
        // cout<<p1.size()<<endl;

	}
	c=q;
	p1 =parent1, p2=parent2;
	// cout<<p1.size()<<endl;
	child2.push_back(c);
	while(p1.size()!=1){
		vector<int>::iterator it1;
		it1 = p1.begin();
		while(*it1 != c){
			it1++;
		}
		vector<int>::iterator it2;
		it2 = p2.begin();
		while(*it2 != c){
			it2++;
		}
		int dx, dy;
		if(it1!=p1.begin()){
			dx = *(it1-1);
		}
		else{
			dx = p1[p1.size()-1]; 
		}
		if(it2!=p2.begin()){
			dy = *(it2-1);
		}
		else{
			dy = p2[p1.size()-1]; 
		}
	
        c = tsp->graph[c][dx] > tsp->graph[c][dy] ? dy:dx;
        p1.erase(it1);
		p2.erase(it2);
        child2.push_back(c);

	}
	

}


void Genetic::doNewCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2){
	int startPos, endPos;
	// do{
		startPos = rand()%parent1.size();
		endPos = rand()%parent1.size();
	// }while(startPos==endPos);
	child1.resize(parent1.size());
	child2.resize(parent2.size());
	for(int i=0;i<child1.size();i++){
		child1[i]=-1;
		child2[i]=-1;
	}
	map<int, int> m;
	for(int i=0;i<child1.size();i++){
		if(startPos<endPos && i>startPos && i<endPos){
			child1[i]=parent1[i];
			m[parent1[i]]=1;
		}
		else if(startPos>endPos){
			if(!(i<startPos && i>endPos)){
				child1[i]=parent1[i];
				m[parent1[i]]=1;
			}
		}
	}

	for(int i=0;i<parent2.size();i++){
		if(m[parent2[i]]!=1){
			for(int j=0;j<child1.size();j++){
				if(child1[j]==-1){
					child1[j]=parent2[i];
					break;
				}
			}
		}
	}

	map<int, int> m1;
	for(int i=0;i<child2.size();i++){
		if(startPos<endPos && i>startPos && i<endPos){
			child2[i]=parent2[i];
			m1[parent2[i]]=1;
		}
		else if(startPos>endPos){
			if(!(i<startPos && i>endPos)){
				child2[i]=parent2[i];
				m1[parent2[i]]=1;
			}
		}
	}

	for(int i=0;i<parent1.size();i++){
		if(m1[parent1[i]]!=1){
			for(int j=0;j<child2.size();j++){
				if(child2[j]==-1){
					child2[j]=parent1[i];
					break;
				}
			}
		}
	}
}

void Genetic::doMCycleCrossover(vector<int>& child1, vector<int>& child2, vector<int> parent1, vector<int> parent2){
	vector<bool> temp1(parent1.size());
	for(int i=0;i<parent1.size();i++){
		temp1[i]=false;
	}
	// cout<<temp1.size()<<endl;
	cout<<"Hi"<<endl;
	child1.resize(parent1.size());
	child2.resize(parent1.size());
	// cout<<child1.size()<<endl;
	int j=0, k=0;
	while(j!=child1.size()){
	// while(parent1.size()!=0){
		int temp=parent2[0];
		int y = -1;
		// cout<<"Hi1"<<endl;
		// cout<<parent1.size()<<" "<<parent2.size()<<endl;
		while(y!=parent1[0] && j!=child1.size()){
			int x = temp;
			child1[j] = temp;
			temp1[temp]=true;
			j++;
			// cout<<j<<endl;
			int i=0;
			for(i=0;i<parent1.size();i++){
				if(x==parent1[i])
					break;
			}
			// cout<<"hi2"<<endl;
			x = parent2[i];
			for(i=0;i<parent1.size();i++){
				if(x==parent1[i])
					break;
			}
			// cout<<"hi3"<<endl;
			child2[k]=parent2[i];
			// cout<<parent2[i]<<endl;
			temp1[parent2[i]]=true;
			k++;
			y = parent2[i];
			x = parent2[i];
			// cout<<"hi4"<<endl;
			for(i=0;i<parent1.size() && j!=child1.size();i++){
				// cout<<"hi5 "<<i<<endl;
				if(x==parent1[i])
					break;
			}
			temp = parent2[i];
			// cout<<"hi6"<<endl;
		}
		if(j!=280){
		vector<int>::iterator it;
		for(it=parent1.begin();it!=parent1.end();){
			if(temp1[*it]){
				vector<int>::iterator it1 = it;
				it++;
				parent1.erase(it1);
			}
			else{
				it++;
			}
		}
		for(it=parent2.begin();it!=parent2.end();){
			if(temp1[*it]){
				vector<int>::iterator it1 = it;
				it++;
				parent2.erase(it1);
			}
			else{
				it++;
			}
		}
		}
	
	}
	cout<<"hey"<<endl;
	cout<<isValidSolution(child1)<<" "<<isValidSolution(child2)<<endl;
}

void Genetic::mutation(vector<int>& child){
	int mutation = rand() % 100 + 1; // random number in [1,100]
	if(mutation <= mutation_rate) // checks if the random number <= mutation rate
	{
		// makes a mutation: change of two genes
		
		int index_gene1, index_gene2;
		index_gene1 = rand() % (tsp->n - 1) + 1;
		index_gene2 = rand() % (tsp->n - 1) + 1;
		
		int aux = child[index_gene1];
		child[index_gene1] = child[index_gene2];
		child[index_gene2] = aux;
		
	}
}

void Genetic::crossOver(vector<int>& parent1, vector<int>& parent2, vector< pair<vector<int> ,int> >& newPopulation)
{
	vector<int> child1, child2;
	
	// map of genes, checks if already are selected
	// doNearestNeighborCrossover(child1, child2, parent1, parent2);
	// doNewCrossover(child1, child2, parent1, parent2);
	// doCrossover(child1, child2, parent1, parent2);
	doMCycleCrossover(child1, child2, parent1, parent2);
	// cout<<"Select the crossover technique:"<<endl;
	// int y;
	// cin>>y;


	mutation(child1);
	mutation(child2);
	
	int total_cost_child1 = isValidSolution(child1);
	int total_cost_child2 = isValidSolution(child2);
	
	// checks if is a valid solution and not exists in the population
	if(total_cost_child1 != -1 && !existsChromosome(child1))
	{
		// add child in the population
		insertBinarySearch(child1, total_cost_child1, newPopulation); // uses binary search to insert
		real_size_population++; // increments the real_size_population
	}
	
	// checks again...
	if(total_cost_child2 != -1 && !existsChromosome(child2))
	{
		// add child in the population
		insertBinarySearch(child2, total_cost_child2, newPopulation); // uses binary search to insert
		real_size_population++; // increments the real_size_population
	}
}


// runs the genetic algorithm
void Genetic::run()
{
	initialPopulation(); // gets initial population
	
	if(real_size_population == 0)
		return;

	for(int i = 0; i < generations; i++)
	{
		vector< pair<vector<int> ,int> > newPopulation;
		vector< pair<int,int> > p;
		newPopulation.push_back(population[0]);

		for(int j=0;j<size_population-1;j++){
			int parent1=0;
			int parent2=0;
			do
			{
				// select two random parents
				parent1 = rand() % size_population;
				parent2 = rand() % size_population;
			}while(parent1 == parent2);
			if(find(p.begin(), p.end(), make_pair(parent1,parent2))==p.end() || j==0){
				p.push_back(make_pair(parent1, parent2));
				crossOver(population[parent1].first, population[parent2].first, newPopulation);
			}
			else
				j--;

		}

		population = newPopulation;
		if(population[0].second < bestPop.second){
			bestPop = population[0];
		}
	}
	
	if(show_population == true)
		showPopulation(); // shows the population
	
	cout << "Best solution: ";
	const vector<int>& vec = bestPop.first;
	// for(int i = 0; i < tsp->n; i++)
	// 	cout << vec[i] << " ";
	// cout << tsp->circuit[0];
	cout << " | Cost: " << bestPop.second;
	string outName = in+".tour";
	ofstream outputStream(outName.c_str());
	// outputStream.open(outFname.c_str(), ios::out);
	outputStream << bestPop.second << endl;
	for (int i=0;i<tsp->n;i++) {
	//for (vector<int>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
		outputStream << vec[i] << endl;
	}
	outputStream<<bestPop.first[0]<<endl;
	//outputStream << *(circuit.end()-1);
	outputStream.close();
}


int Genetic::getCostBestSolution()
{
	if(real_size_population > 0)
		return bestPop.second;
	return -1;
}




//================================ PRINT FUNCTIONS ================================//

void TSP::printResult(){
	string outName = inFile + ".tour";
	ofstream outputStream(outName.c_str());
	// outputStream.open(outFname.c_str(), ios::out);
	outputStream << pathLength << endl;
	for (vector<int>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
	//for (vector<int>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
		outputStream << *it << endl;
	}
	//outputStream << *(circuit.end()-1);
	outputStream.close();
};

void TSP::printPath(){
	cout << endl;
	for (vector<int>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
		cout << *it << " to " << *(it+1) << " ";
		cout << graph[*it][*(it+1)] << endl;
	}
	cout << *(circuit.end()-1) << " to " << circuit.front();

	cout << "\nLength: " << pathLength << endl << endl;
};

void TSP::printEuler() {
	for (vector<int>::iterator it = circuit.begin(); it != circuit.end(); ++it)
		cout << *it << endl;
}

void TSP::printAdjList() {
	for (int i = 0; i < n; i++) {
		cout << i << ": "; //print which vertex's edge list follows
		for (int j = 0; j < (int)adjlist[i].size(); j++) {
			cout << adjlist[i][j] << " "; //print each item in edge list
		}
		cout << endl;
	}
};

void TSP::printCities(){
	cout << endl;
	int i = 0;
	for (vector<City>::iterator it = cities.begin(); it != cities.end(); ++it)
		cout << i++ << ":  " << it->x << " " << it->y << endl;
}