//******************************************************************************************
//
//	Program:	independant_set.cc
//
//	Name:		Max Rego
//	Email: 		mr255509@ohio.edu
//	Date:		9/8/2014
//
//	Project:	HW 1 : Independant Set
//	
//	Comments:	This program will find an Independant Set from a graph
//
//	How to Run: g++ -O2 independant_set.cc
//				./a.out example.dat
//
//********************************************************************************************

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <climits>
#include <stdlib.h> 
#include <algorithm>
using namespace std;

//Prints the Adj List to the screen
void print_AdjList(const vector<set<int> > &adj_list){
	
	for(int vert = 0; vert < adj_list.size(); vert++){  
		for(set<int>::iterator iter=adj_list[vert].begin(); iter!=adj_list[vert].end();++iter){
			cout << *iter << " ";
		}
		cout << endl;
	}
}

//Finds if an element is in the visted/negative set
//input : is the vector negative set, and an int v which is the desired vertex
//output : is 1 if the element is in the negative set or 0 if it is not 
bool find_ele(const vector<int> &negative_set, int v){
	
	if(negative_set.empty()){
		return 0;
	}
	
	vector<int>::const_iterator it = find(negative_set.begin(), negative_set.end(), v);
	if (it != negative_set.end()){
		return 1;
	} else {
		return 0;
	}
}
	
//Returns the vertex with the highest degree that is not in the negative set, may be ties it will choose at least one of them
//Input : Adj List and the Negative Set of visted vertices
//Output : interger of the vertex with at least the highest degree, or -1 if all possibilties are in the negative set ie all has been visited	
int get_high_degree(const vector<set<int> > &adj_list, const vector<int> &negative_set){
	
	int vertex = -1;
	int largest = 0;
	
	for(int i = 0; i < adj_list.size(); i++){
		if (find_ele(negative_set, i) == 0){
			if (adj_list[i].size() >= largest){
				largest = adj_list[i].size();
				vertex = i;
			}
		}
	}
	
	if (vertex > -1){
		return vertex;
	} else {
		return -1;
	}
}

//Returns the vertex with the lowest degree that is not in the negative set, may be ties, it will choose one of them
//Input : Adj List and the Negative Set of visted vertices
//Output : interger of the vertex with at least the highest degree, or -1 if all possibilties are in the negative set ie all has been visited
int get_low_degree(const vector<set<int> > &adj_list, const vector<int> &negative_set){
	
	int vertex = -1;
	int smallest = 1000000; //hard coded to make this easier, should be primed to max value for the graph
	
	for(int i = 0; i < adj_list.size(); i++){
		if (find_ele(negative_set, i) == 0){
			if (adj_list[i].size() <= smallest){
				smallest = adj_list[i].size();
				vertex = i;
			}
		}
	}
	
	if (vertex > -1){
		return vertex;
	} else {
		return -1;
	}
}

//Returns a set with the neighbors to a vertex K
//Inputs : Adj List, Negative set with the visted vertexs, and an interger K with the vertex to find neighbors of
//Output : a Set<int> with all neighbors of a vertex K
set<int> get_linked(const vector<set<int> > &adj_list, const vector<int> &negative_set, int k){
	
	set<int> linked;
	
	for(set<int>::iterator iter=adj_list[k].begin(); iter!=adj_list[k].end();++iter){
		if (find_ele(negative_set, *iter) == 0){
			linked.insert(*iter);
		}
	}
	
	return linked;
}

//Returns a vector with a set of linked neighbors added to a negative set
//Input : Negative_Set with all visted vertex, Linked a set of neighbors to add to the negative set
//Output: A Vector with the negative_set and the linked neighbors combined
vector<int> add_set(const vector<int> &negative_set, const set<int> linked){

	vector<int> linked_added = negative_set;
	for(set<int>::const_iterator i=linked.begin(); i!=linked.end();++i){
		linked_added.push_back(*i);
	}
	return linked_added;
}

//Finds all degree zero vertexs and adds them to my_set and returns the number of them it found
//Input : Adj List, negative set with all visted vertex, and my_set with the final IS
//Output, interger with the number of zero vertex, also adds to my_set for final IS
int get_degree_zero(const vector<set<int> > &adj_list, vector<int> &negative_set, vector<int> &my_set){
	
	int total = 0;
	for(int i = 0; i < adj_list.size(); i++){
		if(find_ele(negative_set, i) == 0){
			set<int> linked = get_linked(adj_list, negative_set, i);
			if (linked.empty()){
				total++;
				my_set.push_back(i);
				negative_set.push_back(i);
			}
		}
	}
	return total;
}

//Finds a single degree 1 vertex, adds it to my_set and adds its neighbors to the negative set
//Input : Adj List, negative set with all visted vertexs, my_set with the final IS
//Output : returns interger 1 if there is a degree 1 vertex not in the negative set, adds that vertex to negative set and my_set
int get_degree_one(const vector<set<int> > &adj_list, vector<int> &negative_set, vector<int> &my_set){
	
	for(int i = 0; i < adj_list.size(); i++){
		if(find_ele(negative_set, i) == 0){
			set<int> linked = get_linked(adj_list, negative_set, i);
			if(linked.size() == 1){
				my_set.push_back(i);
				negative_set.push_back(i);
				set<int>::iterator k = linked.begin(); 
				negative_set.push_back(*k);
				return 1;
			}
		}
	}
	return 0;
}

//Finds a degree two vertex and adds it to my_set and adds its neighbors to the negative set
//Input : Adj List, negative set with all visted vertexs, my_set with the final IS
//returns 1 if there is a vertex of degree 2 in the graph, adds choosen vertex to the negative set and my_set
int get_degree_two(const vector<set<int> > &adj_list, vector<int> &negative_set, vector<int> &my_set){
	
	for(int i = 0; i < adj_list.size(); i++){
		if(find_ele(negative_set, i) == 0){
			set<int> linked = get_linked(adj_list, negative_set, i);
			if(linked.size() == 2){
				my_set.push_back(i);
				negative_set.push_back(i);
				for(set<int>::iterator k=linked.begin(); k!=linked.end();++k){
					negative_set.push_back(*k);
				}
				return 1;
			}
		}
	}
	return 0;
}
	
//Returns the size of largest Independant Set, my_set contains the IS
//Inputs: Adj List passed by reference, negative set to hold the visted vertexs, and my_set to hold the final IS 
int get_IS(const vector<set<int> > &adj_list, vector<int> negative_set, vector<int> &my_set){
	
	//Holds the largest degree vertex remaining in the graph
	int max_vertex = get_high_degree(adj_list, negative_set);
	set<int> linked = get_linked(adj_list, negative_set, max_vertex);
	int max_degree = linked.size();
	
	//Chooses if we are at the base case or if there are still vertexs with degree 3+
	if(max_degree < 3){
		
		int total = 0;
	
		//Remove all degree 1 and 0 vertexs
		while(get_degree_one(adj_list, negative_set, my_set) > 0){
			total++;
		}
		total += get_degree_zero(adj_list, negative_set, my_set);
		
		//find a cycle and remove a vertex, then treat the graph again for snakes and zeros, then checks for another cycle until all are in negative set
		//At this point we can assume that there are only cycles left in the graph if negative set is smaller than adj_list size
		while(negative_set.size() != adj_list.size()){
			int j = get_degree_two(adj_list, negative_set, my_set);
			total++;
			while(get_degree_one(adj_list, negative_set, my_set) > 0){
				total++;
			}
			total += get_degree_zero(adj_list, negative_set, my_set);
		}
		return total;
			
	} else {
		vector<int> negative_set_removed = negative_set;
		negative_set_removed.push_back(max_vertex);
		vector<int> negative_set_choose = add_set(negative_set_removed, linked);
		
		vector<int> my_set_removed = my_set;
		vector<int> my_set_choose = my_set;
		my_set_choose.push_back(max_vertex);
		
		int max_remove_vertex = get_IS(adj_list, negative_set_removed, my_set_removed);
		int max_choose_vertex = get_IS(adj_list, negative_set_choose, my_set_choose);
	
		//choose wheather or not choosing the large vertex or not gives a larger IS
		if(max_remove_vertex > max_choose_vertex){
			
			my_set = my_set_removed;
			return max_remove_vertex;
			
		} else {
			
			my_set = my_set_choose;
			return max_choose_vertex + 1;
				
		}	
	}
}

//Function to confirm the IS is real
//Input : Adju List and the final IS
void confirm_IS(const vector<set<int> > &adj_list, const vector<int> &my_set){
	
	for(int i = 0; i < my_set.size(); i++){
		for(set<int>::iterator k = adj_list[my_set[i]].begin(); k != adj_list[my_set[i]].end(); ++k){
			for(int j = 0; j < my_set.size(); j++){
				if(my_set[j] == *k){
					cout << "Not and Independant Set!!"<<endl<<"********************************************"<<endl;
					return;
				}
			}
		}
	}
	cout << "Confirmed Independant Set!!"<<endl<<"**********************************************"<<endl;
	return;
}

//Main Function
int main(int argc, char *argv[]) {

	//Input Code
	ifstream in1;
	ifstream in2;
	in2.open(argv[1]);
	int m;
	in2 >> m;
	string line;
	getline(in2,line); // Get the newline
	vector<set<int> > adj_list;
	adj_list.resize(m);
	for (int i=0;i<m;i++) {
		getline(in2,line);
		if (!in2.fail()) {
			istringstream in3(line);
			int j;
			in3>> j;
			char c;
			in3 >>c;
			while (!in3.eof()) {
				int k;
				in3 >> k;
				if (!in3.fail()) {
					adj_list[j].insert(k);
				}
			}
		}
	}
	cout << endl<< "**********************************************" <<endl;
	cout << "Number of Vertices = " << adj_list.size() << endl;

	//make containers for negative set and final IS
	vector<int> negative_set;
	vector<int> my_set;
	
	//Find IS and get its size
	cout << "Max Independant Set Size = "<< get_IS(adj_list, negative_set, my_set) << endl << "Max IS = ";

	//Sort the Indepentant Set for printing
	sort(my_set.begin(), my_set.end());

	//print Independant Set
	for(int i = 0; i < my_set.size(); i++){
		cout << my_set[i] << " ";
	}
	cout << endl;
	
	//Confirm IS
	confirm_IS(adj_list, my_set);

}
  
  
  
