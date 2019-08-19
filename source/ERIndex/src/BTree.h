/*
 * BTree.h
 *
 *  Created on: 25/set/2014
 *      Author: fernando
 */

#ifndef BTREE_H_
#define BTREE_H_

namespace std {


// A BTree node
class BTreeNode

{
    private:
        //The simplest B-tree occurs when t = 2 (for t=1 it would be a binary tree)
        int t;

        //Every node can contain at most 2t - 1 keys. Therefore, an internal node can have at most 2t children
        int *keys;

        //Vector of pointers to the children
        BTreeNode **C;

        //Number of keys currently stored in this node
        int n;

        //True if this node is a leaf
        bool leaf;

    public:

        BTreeNode(int t1, bool leaf1);

        // traverse all nodes in a subtree rooted with this node

        void traverse();

        void insertNonFull(int k);

        void splitChild(int i, BTreeNode *y);

        BTreeNode *search(int k);

        friend class BTree;

};


// Class BTree
class BTree
{
    private:
        BTreeNode *root;
        int t;
    public:
        BTree(int _t);
        void traverse();
        BTreeNode* search(int k);
        void insert(int k);
        virtual ~BTree();

};



} /* namespace std */



#endif /* BTREE_H_ */
