/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree search functions

=============================================================================*/


#ifndef SPIDIR_TOP_CHANGE_H
#define SPIDIR_TOP_CHANGE_H

namespace spidir {


void performNni(Tree *tree, Node *nodea, Node *nodeb);
void proposeRandomNni(Tree *tree, Node **a, Node **b);
void performSpr(Tree *tree, Node *subtree, Node *newpos);
void proposeRandomSpr(Tree *tree, Node **subtree, Node **newpos);
bool validSpr(Tree *tree, const Node *subtree, const Node *newpos);


} // namespace spidir

#endif // SPIDIR_TOP_CHANGE_H
