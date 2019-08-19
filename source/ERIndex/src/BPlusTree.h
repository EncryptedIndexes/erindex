/*
 * BPlusTree.h
 *
 *  Created on: 25/set/2014
 *      Author: fernando
 */

#ifndef BPLUSTREE_H_
#define BPLUSTREE_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include "FileBitWriter.h"
#include "BitWriter.h"
#include "FileBitReader.h"
#include <boost/dynamic_bitset.hpp>
#include "EncryptedFileBitWriter.h"
#include "EncryptedFileBitReader.h"

namespace std {

//used to add values and to hold them in memory after the leaf has been deciphered.
//the rangeQuery will process values in this form and it will return them grouped by individualId
typedef struct Value{
	int individualId;
	int factorId;

	Value(){}
	Value(int individualId,int factorId){
		this->individualId=individualId;
		this->factorId=factorId;
	}

	friend bool operator==(const Value& ls, const Value& rs) {
	    return ls.individualId == rs.individualId && ls.factorId == rs.factorId;
	}

	friend ostream& operator<<(std::ostream& out, const Value& v){
	   return out << v.individualId << '/' << v.factorId;
	}

} Value;







//used to store values. A leaf will contain:
// 1) an array of individualsId, encrypted with the system key coupled with iv=nodeNumber
// 2) each value will be stored as a couple <individualRef,factorId>, where individualReference is the
//    position of the individual in the aforementioned array
//
//This strategy produce only a little overhead:
// - the individualsIds will be encrypted in one pass;
// - the referenceIds will be very small integer numbers and they can be encoded with a few number of bits.

//This is acceptable regards to security, because the same individual will be encrypted differently
//in different leafs and so this strategy doesn't reveal to an attacker that all the factors of an individual (even not knowing who he is)
//start in certain reference positions (posTree) or that they can be reconstructed from certain reverse suffix array indexes
//(reverseTree) or ...
typedef struct StoredValue{
	int individualRef;
	int factorId;
} StoredValue;




class BPlusTree;

/* the tree directory must be stored at the end of the bytestream because
 * we know the compressed starts positions only after the compression of each node
 */
typedef struct TreeDirectory{
	//stored informations
	uint treeOffset; //offset of the tree in the LZ-Index bytestream
	        		 //(this must be used to obtain the absolute offset of the tree nodes and directory
	                 //within the whole LZ-Index)
	uint nodesNumber;
	uint lastLeafNodeNumber;

	//compressed nodes starts positions bitmap

	//runtime-only informations
	vector<uint64_t> nodesStartsPositions;  //start positions of all the nodes following the first
	boost::dynamic_bitset<> nodesStartsPositionsBitmap;  //nodes starts positions in the bytestream
	boost::dynamic_bitset<> nodesStartsPositionsCompressedBitmap;
} TreeDirectory;

// A B+ Tree node
class BPlusTreeNode{
    private:
        //The simplest B-tree occurs when t = 2 (for t=1 it would be a binary tree)
        int t;
        //Every node can contain at most 2t - 1 keys. Therefore, an internal node can have at most 2t children
        int *keys;
        // Values corresponding to the above keys;
        vector<Value> *values; //each element of this array is the array of values corresponding to the key in the same position
        //Vector of pointers to the children
        BPlusTreeNode **C;
        //Vector of disk pointers (node numbers) to the children
        int *CNN;  //children node numbers
        //Number of keys currently stored in this node
        int n;
        //True if this node is a leaf
        bool leaf;
        BPlusTreeNode* parent;
        BPlusTreeNode* previousSibling;
        BPlusTreeNode* nextSibling;
        BPlusTree* tree;
        uint nodeNumber;

        //association between individualId and individualRef, used while building the index:
        // first: individualId, second:individualReference
        map<int,int> individualRefs;
        map<int,EncryptionContext *> individualEncryptionContexts;

        bool binarySearchInNode(int k, int *position);
        void saveKeys(FileBitWriter &bitWriter);
        void loadKeys(FileBitReader *bitReader);
    public:
        BPlusTreeNode(int t1, bool leaf1);
        void nodeStats(int *leafNodes, int *indexNodes);
        // traverse all nodes in a subtree rooted with this node
        void traverse();
        void insert(int k,BPlusTreeNode *c);
        void insert(int k,Value value);
        void split();
        bool search(int k,vector<Value> *value);
        BPlusTreeNode *searchInsertionLeaf(int k);
        void dumpSubTree(int level);
        void save(EncryptedFileBitWriter &bitWriter);
        void load(EncryptedFileBitReader *bitReader);
        friend class BPlusTree;
};

class LZIndex;

// Class BPlusTree (B+ tree)
class BPlusTree{
    private:
		LZIndex *index;
		uint nodesNumber;
		uint savedNodes;
        BPlusTreeNode *root;
        int t;
        BPlusTreeNode *firstLeaf;
        BPlusTreeNode *lastLeaf;

        int referenceLength; //length of the reference genome
        int numberOfIndividuals; //number of individuals in the LZ-index
        int maximumIndividualId;
        int maximumNumberOfFactorsPerIndividual; //maximum number of factors, computed over all the individual factorizations
        uint64_t treeDirectoryOffset; //offset of the tree directory in the BPlusTree bytestream
        TreeDirectory  treeDirectory; //treeDirectory
        EncryptedFileBitReader *bitReader;
        uint32_t baseNonce;  //number used as base for all the nonces used to encrypt this tree

        //TODO: replace the following vector with a LRU cache
        vector<BPlusTreeNode*> nodes;

        BPlusTreeNode* searchInsertionLeaf(int k);
    public:
        BPlusTree(LZIndex *index,int _t,uint64_t referenceLength);
        BPlusTreeNode *newNode(bool leaf);
        void traverse();
        void treeStats(int* leafNodes,int *indexNodes);
        void treeTraverse();
        void insert(int k,Value value);
        void dumpTree();
        virtual ~BPlusTree();
        bool search(int k,vector<Value>* value);
        bool isValueInRange(int l,int u,Value value);
        void load(EncryptedFileBitReader *bitReader,uint64_t treeDirectoryOffset);
        void loadAllNodesInMemory();
        uint64_t save(EncryptedFileBitWriter &bitWriter,uint64_t *treeDirectoryOffset);
        void setBaseNonce(uint32_t baseNonce);
        void setNumberOfIndividuals(int numberOfIndividuals);
        void setMaximumNumberOfFactorsPerIndividual(int maximumNumberOfFactorsPerIndividual);
        void setMaximumIndividualId(int maximumIndividualId);
        map<int,vector<pair<int,int>>> rangeQuery(int l,int u);
        vector<pair<int,vector<Value>>> rangeQueryAllInMemory(int l,int u);
        BPlusTreeNode *getNode(int nodeNumber);
        friend class BPlusTreeNode;
};


} /* namespace std */

#endif /* BPLUSTREE_H_ */
