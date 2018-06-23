#ifndef TREE_H
#define TREE_H
#include "box.h"

template <class T>
class Tree
{
public:
	Tree(void);
	Tree(const T &data);
	Tree(const T &data, Tree<T> *left, Tree<T> *right);
	
	T data;
	Tree<T>* left;
	Tree<T>* right;

	void split(const int &minTriangles, const Mesh &mesh);

	void print(const int depth);

	std::vector<Box>* getLeaves(void);

	void Tree<Box>::highlightEdges(void);
};


#endif