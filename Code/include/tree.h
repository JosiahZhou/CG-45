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

	int size(void)
	{
		if (this == NULL) return 0;
		return this->left->size() + this->right->size() + 1;
	}

	void split(const int minTriangles, const int maxLevel, Mesh &mesh);

	void print(const int depth);

	std::vector<Box*> getLeaves(void);

	void highlightEdges(void);
};
#endif