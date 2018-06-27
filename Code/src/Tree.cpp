#include "tree.h"


template<class T> Tree<T>::Tree(void)
{
	this->data = 0;
	this->left = 0;
	this->right = 0;
}

template<class T> Tree<T>::Tree(const T &data)
{
	this->data = data;
	this->left = 0;
	this->right = 0;
}

template<class T> Tree<T>::Tree(const T &data, Tree<T>* left, Tree<T>* right)
{
	this->data = data;
	this->left = left;
	this->right = right;
}

// Splits the box of the tree into two and allocates them to the left
// and right trees until the minimum triangles in a box has been
// exceeded.
void Tree<Box>::split(const int &minTriangles, const Mesh &mesh)
{
	if (data.triangles.size() < minTriangles) return;
	std::pair<Box, Box> boxes = data.split(mesh);

	this->left = new Tree<Box>(boxes.first);
	this->right = new Tree<Box>(boxes.second);

	this->left->split(minTriangles, mesh);
	this->right->split(minTriangles, mesh);
}

// https://stackoverflow.com/questions/13484943/print-a-binary-tree-in-a-pretty-way
int rec[1000006];
// Prints the BoxTree in directory-format
// TODO: fix this
void Tree<Box>::print(const int depth)
{
	int i;
	if (this == NULL)return;
	printf("\t");
	for (i = 0; i < depth; i++)
		if (i == depth - 1)
			printf("%s---", rec[depth - 1] ? "|" : "|");
		else
			printf("%s   ", rec[i] ? "|" : "  ");
	printf("%ld\n", this->data.triangles.size());
	rec[depth] = 1;
	this->right->print(depth + 1);
	rec[depth] = 0;
	this->left->print(depth + 1);
}

// Getter for all the leaf boxes inside the Tree<Box>. Boxes are added recursively to a tree.
std::vector<Box>* Tree<Box>::getLeaves(void)
{
	if (this == NULL) return &std::vector<Box>();
	if (this->left == NULL && this->right == NULL) return &std::vector<Box>{ this->data };
	std::vector<Box>* boxes = this->left->getLeaves();
	std::vector<Box>* rightBoxes = this->right->getLeaves();
	boxes->insert(boxes->end(), rightBoxes->begin(), rightBoxes->end());
	return boxes;
}

// Highlight all the edges of the boxes in the tree.
void Tree<Box>::highlightEdges(void)
{
	this->data.highlightEdges();
	if (left == NULL) return;
	this->left->highlightEdges();
	this->right->highlightEdges();
}