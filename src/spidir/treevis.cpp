/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Visualize trees

=============================================================================*/

// c++ headers
#include <assert.h>
#include <stdio.h>

// spidir headers
#include "Tree.h"
#include "common.h"
#include "Matrix.h"


namespace spidir {


void drawLine(Matrix<char> &matrix, char chr, int x1, int y1, int x2, int y2)
{
    float stepx, stepy;
    int steps;
    
    steps = max(abs(x2-x1), abs(y2-y1));
    
    stepx = float(x2 - x1) / steps;
    stepy = float(y2 - y1) / steps;
    
    float x = x1;
    float y = y1;
    for (int i=0; i<steps+1; i++) {
        matrix[int(y)][int(x)] = chr;
        x += stepx;
        y += stepy;
    }
}


void drawText(Matrix<char> &matrix, const char *text, int x, int y)
{
    int len = strlen(text);
    
    for (int i=0; i<len; i++) {
        matrix[y][x+i] = text[i];
    }
}


// TODO: finish display
void displayNode(Matrix<char> &matrix, Node *node, int *xpt, int* ypt)
{
    int x = xpt[node->name];
    int y = ypt[node->name];
    const int leadSpacing = 2;

    
    for (int i=0; i<node->nchildren; i++) {
        displayNode(matrix, node->children[i], xpt, ypt);
    }

    // horizontal line
    if (node->parent) {
        drawLine(matrix, '-', x, y, 
                 xpt[node->parent->name], y);
    }

    // vertical line
    if (!node->isLeaf()) {
        int l = node->nchildren - 1;
        drawLine(matrix, '|', x,
                              ypt[node->children[0]->name],
                              x,
                              ypt[node->children[l]->name]);
        matrix[ypt[node->children[0]->name]][x] = '/';
        matrix[ypt[node->children[l]->name]][x] = '\\';

        matrix[y][x] = '+';                    
    } else {
        drawText(matrix, node->longname.c_str(), x + leadSpacing, y);
    }
}


int treeLayout(Node *node, int *xpt, int *ypt,
               float xscale, int yscale, int y=0)
{
    if (node->parent == NULL) {
        xpt[node->name] = int(xscale * node->dist);
        ypt[node->name] = 0;
    } else {
        xpt[node->name] = xpt[node->parent->name] + int(xscale * node->dist + 1);
        ypt[node->name] = y;
    }
    
    if (node->isLeaf()) {
        y += yscale;
    } else {
        for (int i=0; i<node->nchildren; i++) {
            y = treeLayout(node->children[i], xpt, ypt, xscale, yscale, y);
        }
        int l = node->children[node->nchildren-1]->name;
        ypt[node->name] = (ypt[l] + ypt[node->children[0]->name]) / 2;
    }
    
    return y;
}


// find maximum width of layout
int treeLayoutWidth(Tree *tree, int *xpt, int *ypt, int labelSpacing=2)
{
    int width = 0;    
    for (int i=0; i<tree->nnodes; i++) {
        int extra = tree->nodes[i]->isLeaf() ?
            labelSpacing + tree->nodes[i]->longname.size() : 0;
        if (xpt[i] + extra > width) 
            width = xpt[i] + extra;
    }
    return width;
}


// TODO: add names
void displayTree(Tree *tree, FILE *outfile, float xscale, int yscale)
{
    const int labelSpacing = 2;

    ExtendArray<int> xpt(tree->nnodes);
    ExtendArray<int> ypt(tree->nnodes);
    treeLayout(tree->root, xpt, ypt, xscale, yscale);
    int width = treeLayoutWidth(tree, xpt, ypt, labelSpacing);
    
    // init character matrix with blanks
    Matrix<char> matrix(tree->nnodes+1, width+1);
    matrix.setAll(' ');
    
    displayNode(matrix, tree->root, xpt, ypt);
    
    // write out matrix
    for (int i=0; i<matrix.numRows(); i++) {
        for (int j=0; j<matrix.numCols(); j++)
            fprintf(outfile, "%c", matrix[i][j]);
        fprintf(outfile, "\n");
    }
}

void displayTreeMatrix(Tree *tree, float xscale, int yscale, 
                       char ***matrix, int *nrows, int *ncols)
{
    const int labelSpacing = 2;

    ExtendArray<int> xpt(tree->nnodes);
    ExtendArray<int> ypt(tree->nnodes);
    treeLayout(tree->root, xpt, ypt, xscale, yscale);
    int width = treeLayoutWidth(tree, xpt, ypt, labelSpacing);
    
    Matrix<char> matrix2(tree->nnodes+1, width+1);
    matrix2.setAll(' ');
    
    displayNode(matrix2, tree->root, xpt, ypt);
    
    // copy text to new allocated matrix
    *nrows = matrix2.numRows();
    *ncols = matrix2.numCols();
    *matrix = new char* [*nrows];
    
    for (int i=0; i<matrix2.numRows(); i++) {
        (*matrix)[i] = new char [*ncols+1];
        for (int j=0; j<matrix2.numCols(); j++)
            (*matrix)[i][j] = matrix2[i][j];
        (*matrix)[i][*ncols] = '\0';
    }
}


} // namespace spidir
