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

void Tree<Box>::split(const int &minTriangles, const Mesh &mesh)
{
	std::pair<Box, Box> boxes = data.split(minTriangles, mesh);

	this->left = new Tree<Box>(boxes.first);
	this->right = new Tree<Box>(boxes.second);

	if (left == 0)
	{
		this->left->split(minTriangles, mesh);
		this->right->split(minTriangles, mesh);
	}
}

// https://stackoverflow.com/questions/13484943/print-a-binary-tree-in-a-pretty-way
int rec[1000006];
// Prints the BoxTree in directory-format
void Tree<Box>::print(const int depth)
{
	if (this == NULL) return;

	printf("%ld\n", this->data.triangles.size());
	printf("|      ");

	this->right->print(depth + 1);
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