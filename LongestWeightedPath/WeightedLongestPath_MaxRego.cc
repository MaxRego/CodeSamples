//*******************************************************************
// File: WeightedLongestPath_MaxRego.cc
// Author: Max Rego
// Purpose:  This example code reads in a tree with weights on the edges.
// The trees are stored using left-child/right sibling notation.   
// Each tree node has a pointer to the parent of that node as well.
//
// Notes:
//	Possible Longest Paths:
//		1. A Path including a parent and two of its children (upsidedown U Shape)
//		2. A Path including a parent as an endpoint and some node in the tree (straight path)
//	Compile with:
//		g++ read_tree.cc
//	Run With:
//		./a.out < tree.dat
//
//********************************************************************** 
#include <iostream>
#include <cassert>
#include <vector>
#include <stack>
#include <string>
using namespace std;

// Stores a tree using left-child, right sibling notation
class Tree_Node {
public:
	string name;
	Tree_Node *parent;
	Tree_Node *lc;
	Tree_Node *rs;
	int weight;
};

// Stores a solution to the longest path problem
class Solution {
public:
	string left;
	string right;
	int weight;
};

//**********************************************************************
Solution *sol = new Solution; 			//Saves the longest Upsidedown U shape
Solution *solLinePath = new Solution;	//Saves the longest Straight Path
//**********************************************************************

//**********************************************************************
// Function: insert_right
// Purpose: Inserts the node rs into the tree point given by the left-child 
// pointer lc.  If lc is NULL, lc becomes rs, and we are done.  If 
// ls is not NULL, rs becomes the right most node in this chain.
//************************************************************************
void insert_right(Tree_Node *parent, Tree_Node *&lc, Tree_Node *rs) {
	if (lc == NULL) {
		lc = rs;
		return;
	}
  
	Tree_Node *ns = lc;
	while (ns->rs !=NULL) { 
		ns = ns->rs;
	}
	ns->rs = rs;
	rs->parent = parent;
}
//*********************************************************************
// Function: read_tree
// Purpose: Reads a tree from the stream in (cin in this program), and 
// store that tree in the node returned by this function.
//**********************************************************************
Tree_Node *read_tree(Tree_Node *parent, istream &in) {
	char c;
	in >> c;  // Get the next character.
	assert (c == '(');
	string t;
	in >> t;
	int w;
	in >> w;
	//  char c;
	Tree_Node *temp;
	temp = new Tree_Node;
	temp->name = t;
	temp->weight = w;
	temp->lc = NULL;
	temp->rs = NULL;
	temp->parent = parent;
	in >> c;
	//cout << "* " << c << endl;
	//cout << "read " << t << endl;
	while (c == '(') {
		//cout << "Read Clild for " << temp->name << endl;
		in.putback(c);
		Tree_Node *rs; 
		rs = read_tree(temp,in);
		insert_right(temp,temp->lc,rs);
		in >> c;
		//cout << "** " << c << endl;
	}
	if (c == ')') {
		// Good !
		//cout << "Here " << endl;
		//cout << temp->name << endl;
		return temp;
   
	} else {
		cout << "Invalid Tree" << endl;
		return NULL;
	}
}

//**************************************************************************
// Function: Print_Tree
// Purpose: Prints the tree given by root to the standard out (cout)
//**************************************************************************
void Print_Tree(Tree_Node *root) {
	
	cout << "(";
	cout << root->name;
	cout << " " << root->weight;
	Tree_Node *t = root->lc;
	while (t!=NULL) {
		Print_Tree(t);
		t=t->rs;
	}
	cout << ")";
}

//**************************************************************************
// Function: getPath
// Purpose: to find the longest weighted path in the tree in polynomial time
// Comments:
//
//	The interger returned by this funtion is the length of the longest path leaving the root,
//	ie the path that has the root as an endpoint, not in the middle
//
//	The class Solution holds the longest path that is in the subtree or includes the root in
//	the middle of the path
//
//***************************************************************************
int getPath(Tree_Node *root, string &lEnd, string &sLEnd){

	Tree_Node *tmp = root -> lc;
	int lPath = -99999;
	int sLPath = -99999;
	string lName = "NULL";
	string sLName = "NULL";
	
	//Base Case, No children at all
	if (tmp == NULL){
		lEnd = root -> name;
		sLEnd = sLName;
		return root -> weight;

	//Has at least one child	
	} else {
		//get left most child info
		lPath = getPath(tmp, lName, sLName);
		tmp = tmp -> rs;
		//If root has only one child
		if(tmp == NULL){
			//if the longest path is negative dont include it
			if(lPath < 0){
				lEnd = root -> name;
				sLEnd = sLName;
				return root-> weight;
			//if the longest path is positive include it in the path
			} else {
				//check if root connection is negative for straight path inside the tree
				if(root -> weight < 0){
					if(lPath > solLinePath->weight){ 
						solLinePath -> left = lName;
						solLinePath -> right = root->name;
						solLinePath -> weight = lPath;
					}
				}
				lEnd = lName;
				sLEnd = sLName;
				return lPath + root-> weight;
			}
		}
		//has at least two children	
		while (tmp != NULL) {
			//get next childs longest path
			string lNameTmp = "NULL";
			string sLNameTmp = "NULL";
			int tmpPath = getPath(tmp, lNameTmp, sLNameTmp); 
			//compare it to the two longest paths off of the root
			if (tmpPath > lPath){
				sLPath = lPath;
				sLName = lName;
				lPath = tmpPath;
				lName = lNameTmp;
			} else if (tmpPath > sLPath){
				sLPath = tmpPath;
				sLName = lNameTmp;
			}
			tmp = tmp -> rs;
		}
	}
	
	//calculate the longest U shape and compare it to the saved longest
	int includingRoot = lPath + sLPath;
	if(includingRoot > sol-> weight){
		sol -> weight = includingRoot;
		sol -> left = lName;
		sol -> right = sLName;
	}
	
	//if the root connection is negative then save the straight path to this root node
	if(root -> weight < 0){
		//only save it if the longest path below this root is positive
		if(lPath < 0){
			lEnd = root -> name;
			sLEnd = sLName;
			return root -> weight;
		} else {
			//save the longest straight path if this one is longer than the saved one
			if(lPath > solLinePath->weight){ 
				solLinePath -> left = lName;
				solLinePath -> right = root->name;
				solLinePath -> weight = lPath;
			}
		}
	}

	//return the longest path catch-all return
	lEnd = lName;
	sLEnd = sLName;
	return lPath + root -> weight;
	
}

int main() {
	
	//holds the tree
	Tree_Node *root;
	
	//primes the weight of the solutions
	sol -> weight = -9999999;
	solLinePath -> weight = -999999;
	
	//read in the tree
	root = read_tree(NULL,cin);

	//saves the longest straight path including the root of the whole tree
	string lEnd;
	string sLEnd;
	
	//get the longest path
	int k = getPath(root, lEnd, sLEnd);

	//print results
	cout << "*******************************************************************" << endl;
	cout << "Tree : "; 
	Print_Tree(root);
	cout << endl;
	if(k < 0 && sol->weight < 0 && solLinePath -> weight < 0){
		cout << "The Longest Path in the Tree is : Weight: 0, Longest Path is the Root Back to Itself" <<endl;
	} else if(k > sol -> weight){
		if(k > solLinePath->weight){
			cout << "The Longest Path in the Tree is : Weight: " << k << " Endpoints: " << root -> name << " -> "  << lEnd << endl;
		} else if (solLinePath->weight > sol-> weight){
			cout << "The Longest Path in the Tree is : Weight: " << solLinePath -> weight << " Endpoints: " << solLinePath -> left << " " << solLinePath -> right << endl;
		}
	} else if (sol -> weight > solLinePath -> weight){
		cout << "The Longest Path in the Tree is :  Weight: " << sol -> weight << "  Endpoints: " << sol -> left << " " << sol -> right << endl;
	} else {
		cout << "The Longest Path in the Tree is : Weight: " << solLinePath -> weight << " Endpoints: " << solLinePath -> left << " " << solLinePath -> right << endl;
	}
	
	cout << endl << endl;
	cout << "Other Interesting Calculations: " << endl;
	cout << "Longest Straight Path Including the Root of the Whole Tree :   Weight: " << k << " Endpoint: " << root -> name <<" "<< lEnd << endl;
	cout << "Longest Upsidedown U Shape Longest Path (may or may not include root) :  Weight: " << sol -> weight << "  Endpoints: " << sol -> left << " " << sol -> right << endl;
	cout << "Longest Straight Path That May Be Inside the Tree (May Not Include Root of Tree) : Weight: " << solLinePath -> weight << " Endpoints: "; 
	cout << solLinePath -> left << " " << solLinePath -> right << endl;
	cout << "*******************************************************************" << endl;
	
}

