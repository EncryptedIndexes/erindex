/*
 * LZIndex.h
 *
 *  Created on: 21/set/2014
 *      Author: fernando
 */

#ifndef LZINDEX_H_
#define LZINDEX_H_

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <bitset>
#include "LZFactor.h"
#include "BPlusTree.h"
#include "Database.h"
#include "Individual.h"
#include <boost/dynamic_bitset.hpp>
#include <sdsl/int_vector.hpp>
#include "EncryptedFileBitWriter.h"
#include "EncryptedFileBitReader.h"
extern "C"{
	#include "../srcFMIndex/common.h"
}

#define BTREE_HALF_PAGE_SIZE 16

namespace std {


typedef struct Reference{
	string referenceIdentifier;
	string reverseIndexFileName;
	string forwardIndexFileName;
	string correspondenceFilesName;
} Reference;


typedef struct Header{
	//stored informations
	int referenceIdentifier;   //identifier of the reference sequence (32 bit)
	int numberOfIndividuals; //number of individuals in the index
	int *individualIds; //individual Ids in order of factorization
	int maximumFactorsNumber; //maximum number of factors in a factorization
	int blockSize;     //size of each block (number of the factors contained in each block) (32 bit)
	int lmax; //overall factors maximum length (32 bit)
	int reverseTreeDirectoryOffset; //(32 bit)
	int forwardTreeDirectoryOffset; //(32 bit)
	int posTreeDirectoryOffset; //(32 bit)
	int *factorizationStartOffsets;	//array containing the single factorization starts positions within the compressed bytestream

	//here we'll store (over a byte) the number of bits needed to encode the char map
	//here we'll store the char map

	//here we'll store the blockStarts array
	//runtime-only informations
	map<int,int> factorizationIndexes;
	int maximumIndividualId; //maximum individualId (among those contained in the index)
	boost::dynamic_bitset<> charMap; //bitmap used to identify the characters occurring in T
	boost::dynamic_bitset<> compressedCharMap;  //the above charMap, compressed with RLE0
	//int blocksNumber; //number of factor blocks in the index

	uint64_t save(EncryptedFileBitWriter &bitWriter,int *suspendedInfosStartPosition,
			EncryptionContext **suspendedInfosEncryptionContext,int *cRemappings,vector<int> &cInverseRemappings);
	void load(BitReader &bitReader, int *cRemappings,vector<int> &cInverseRemappings);
} Header;








typedef struct Block{
	//stored informations
	//here we'll store the absolute text position of the first factor (over the number of bit required to encode the textLength)
	//here we'll store the size in bits of the compressed relative positions bitmap
	//here we'll store the compressed positions bitmap;
	//here we'll store the saiStarts array

	//runtime-only informations
	/*
	 * L'obiettivo è codificare anche i fattori-extra senza aggiungere un solo bit alla codifica dei fattori standard, che costituiscono di gran lunga la norma (da qui il loro nome).
	   Detto ecn (extra-characters number) il numero dei caratteri extra effettivamente contenuti nel file, ciascuno di essi verrà codificato con un numero [0,ecn-1]
	   che corrisponde al suo codice nell'alfabeto compatto costituito dai soli caratteri extra.
	   Le suffix array positions, che variano nel range [-1,referenceLength[ verranno, quindi, memorizzate con numeri nell'intervallo [ecn,referenceLength+ecn], in modo che ecn corrisponda a -1,ecn+1 a 0
	   e così via. Questa strategia di memorizzazione è molto efficace, in quanto, essendo molto basso il numero degli extra-characters contenuti nel testo, è estremamente improbabile che incrementi il numero di bit necessario a memorizzare le suffix array positions.
	 */

	//LZIndex *index;
	Factorization *factorization;

	int firstFactorAbsolutePosition;
	int relativePositionBits;
	boost::dynamic_bitset<> positionsBitmap;  //bitmap storing the starting positions of the factors in T,
									           //relatively to the starting position of the first factor in this block
	boost::dynamic_bitset<> compressedPositionsBitmap;  //the compressed positions bitmap
	int* saiStarts; //contains, for each standard factor, the suffix array starting position in the backward search
				    //on the reverse index and, for each extra-factor, the extra-character.
	char *letters;

	int factorsNumber;  //number of factors in this block
	vector<LZFactor*> factors;
	int blockIndex; //block index


	Block(Factorization *factorization,int blockIndex);
	int getEstimatedSize(int blockSize,int textLength,int referenceLength);
	//int getFactorStartPosition(int factorOffset);
	LZFactor* getFactor(int factorOffset);

	uint64_t save(BitWriter &bitWriter, int textLength, int referenceLength,int factorsNumber,int blockSize,int compactAlphabetSize);
	void load(BitReader &bitReader, int textLength,  int referenceLength,int totalFactorsNumber,int blockSize,int compactAlphabetSize);
} Block;


typedef struct FactorizationHeader{
	//stored informations
	int textLength;	   //length of the factorized text T (32 bit)
	int factorsNumber;
	int blockStartsOffset;

	//runtime-only informations
	int blocksNumber;//number of factor blocks for the individual
};

class Database;

//An individual sequence factorization: it's stored as an header followed by a sequence of blocks
typedef struct Factorization{
	//stored informations
	FactorizationHeader header;
	vector<Block*> blocks;
	//runtime-only informations
	LZIndex &index;
	Individual *individual;
	vector<LZFactor*> factors;
	int *blockStarts;  //array containing the block starts respect to the first byte of this factorization


	Factorization(LZIndex &index,Individual *individual);
	inline Block* getBlock(int blockIndex,EncryptedFileBitReader &bitReader);
	inline LZFactor* getFactor(int factorIndex,EncryptedFileBitReader &bitReader);
	uint64_t save(EncryptedFileBitWriter &bitWriter);
	void load(EncryptedFileBitReader &bitReader);
	void loadAllBlocks(EncryptedFileBitReader &bitReader);
	friend class Database;
} Factorization;


typedef struct SearchStatistics{
	double locateTime;  				//milliseconds
	double internalLocateTime;
	double overlappingLocateTime;
	double leftOrRightSideLocateTime;
	double assemblingTime;
	bool withExtraCharacters;
} SearchStatistics;



typedef struct ExtraSequence{
	unsigned int position;
	unsigned int length;
	unsigned int lsLength;
	unsigned int rsLength;
	string sequence;
	string ls;
	string rs;
} ExtraSequence;


typedef struct Occurrence{
	int factorIndex;
	int factorOffset;
	int endingFactorIndex;
	int endingFactorOffset;
	int textPosition;

	bool operator < (const Occurrence& other) const
	    {
	        return (factorIndex < other.factorIndex  ||
	        		(factorIndex == other.factorIndex && factorOffset < other.factorOffset));
	    }
} Occurrence;





typedef struct IndividualResult{
	int64_t individualId;
	vector<Occurrence> occurrences;
};




class LZIndex {
public:
	LZIndex(Database &database);
	LZIndex(Database &database,int referenceIdentifier,int blockSize);
	virtual ~LZIndex();
	void open(string indexFileName);
	void save(string indexFileName);
	void dumpIndex(string dumpFileName);
	void close();

	void loadAllInMemory();
	void loadReferenceInMemory();
	void buildIndex(int numberOfIndividuals,Individual **individuals,string *fastaFileNames);
	Factorization *addToIndex(Individual *individual,string fastaFileName,int lookAheadWindowSize,boost::dynamic_bitset<> &threadCharMap,int &threadLMax);

	vector<int64_t> getIndividualIds();
	void doLZFactorization(Factorization *f,const char *buffer,int len,boost::dynamic_bitset<> &threadCharMap, int &threadLMax);
	void locate(const string &pattern,SearchStatistics &stats,map<int64_t,vector<Occurrence*>> &results);

	void findLongestFactor(const char *buffer,int len,int i,int h,int &longest_i,int &longest_l,int &longest_sai_rev,int &longest_numberOfN,
			boost::dynamic_bitset<> &threadCharMap);
	void noLookAheadFactorization(Factorization *f,const char *buffer,int len,int &currentTextPosition,int &currentFactorIndex,
			boost::dynamic_bitset<> &threadCharMap);
	Factorization *getFactorization(int64_t individualId);
	int getFactorsNumber(int64_t factorizationIndex);
	LZFactor* getFactor(int64_t factorizationIndex,int factorIndex);
	string getSequence(int64_t individualId);
	vector<LZFactor*> getFactors(int64_t factorizationIndex);
	string getFactorString(int suffixArrayPosition,int length);
	void setBlockSize(int blockSize);
	int getBlockSize();
	friend Block;
	friend Factorization;
	friend BPlusTree;
	friend BPlusTreeNode;
private:
	void buildSubstringsOccsMap();
	void checkReverseReferenceIndex();  //ensures that the reverse reference index is opened
	void openReverseReferenceIndex();
	void checkForwardReferenceIndex();  //ensures that the reference index is opened
	void checkForwardSuffixArray();
	void openForwardReferenceIndex();
	void closeForwardReferenceIndex();
	void closeReverseReferenceIndex();
	void checkCorrespondenceFiles();  //ensures the r2f and f2r files (establishing the correspondence between the reverse and the forward index) are opened
	void openCorrespondenceFiles();
	void openForwardSuffixArray();

	//Block* getBlock(int blockIndex);
	//bool checkPatternCharacters(const string &pattern);
	void locateInternalOccs(const string &pattern,map<int64_t,vector<Occurrence*>> &results);
	void locateOverlappingOccs(const string &pattern,map<int64_t,vector<Occurrence*>> &results);
	int findLeftSideLongestSuffix(string leftSide);
	int findRightSideLongestPrefix(string rightSide);
	inline map<int,vector<int>> findLeftSideFactors(string ls,int *lslsLength);
	inline map<int,map<int,vector<int>>> findLeftSideFactorsEC(string ls);
	inline map<int,vector<int>> findRightSideFactors(string rs,int *rslpLength);
	inline vector<int> nestedLoopsJoin(vector<int> &leftSideFactors,vector<int> &rightSideFactors,int joinStep);
	bool verifyPatternRemainingPart(Factorization *fn,string pattern,int splitPoint,Occurrence *occurrence,int lsVerifiedLength,int rsVerifiedLength);

	//void loadCorrespondenceArrayFromDisk(FILE *f,int *a,int n);
	void loadCorrespondenceArrayFromDisk(string fileName,int *a,int n);
	void findStartingAndPrefixSuffix(int sai_rev, int l,int *sai_rev_start,int *sai_pref,int *straightTextPosition);
	void forwardBackwardSteps(int sai,int numOfSteps,int *sai_start,int *sai_pref);
	int go_back(FM_INDEX *Infile,int row, int len, char *dest, bwi_out *s);
	void fbit_flush_simple(FILE *f,int* Bit_buffer_size, uint32* Bit_buffer);
	void fbit_write24_simple(FILE *f,int* Bit_buffer_size, uint32* Bit_buffer, int n, int vv);
	void fbit_write_simple(FILE *f, int* Bit_buffer_size, uint32* Bit_buffer,int n, int vv);

	int fbit_read24_simple(FILE *f, int* Bit_buffer_size, uint32* Bit_buffer, int n);
	int fbit_read_simple(FILE *f, int* Bit_buffer_size, uint32* Bit_buffer, int n);
	bool backwardStepSuccessful(FM_INDEX *referenceIndex, bwi_out *s,
			unordered_map<char, unordered_map<int,int>> *occsCache,
			char c,int sp,int ep,int *trySp,int *tryEp);
	map<int,vector<int>> getFactorsInRange(BPlusTree *tree, int sp,int ep);
	bool searchPatternInReferenceIndex(FM_INDEX *referenceIndex, bwi_out *s,
				unordered_map<char, unordered_map<int,int>> *occsCache,
				string pattern,int *sp,int *ep);
	bool searchPatternReverseInReferenceIndex(FM_INDEX *referenceIndex, bwi_out *s,
					unordered_map<char, unordered_map<int,int>> *occsCache,
					string pattern,int *sp,int *ep);
	inline bool checkPatternCharacters(const string &pattern);
	inline bool checkPatternCharacters(const string &pattern,bool *withExtraCharacters,ExtraSequence *extraSequence);
	inline void locateNoExtraCharacters(const string &pattern,map<int64_t,vector<Occurrence*>> &results,SearchStatistics &stats);
	inline void locateExtraCharacters(const string &pattern,ExtraSequence &extraSequence,map<int64_t,vector<Occurrence*>> &results,SearchStatistics &stats);

	string loadFasta(string fastaFileName);
	string loadFasta(string fastaFileName,int size);
	bool tryMatch(LZFactor *factor,int factorOffset,string pattern, int actualPosition,unsigned int &steps);
	bool tryBackwardMatch(LZFactor *factor,int factorOffset, string pattern, int actualPosition,unsigned int &steps);
	inline int reverseBackwardStep(int sai, char *character);

	void factorizationThread(int threadNumber);

	Database &database;
	int referenceIdentifier;
	string indexFilename;
	string forwardReferenceIndexFileName;
	string reverseReferenceIndexFileName;
	string correspondenceFilesName;
	string forwardSuffixArrayFileName;  // usato per rimediare alla scarsa efficienza di FM-INDEX v1 nel marcare le righe
	FM_INDEX *forwardReferenceIndex;
	FM_INDEX *reverseReferenceIndex;
	//FILE *forwardSAFile;
	//int *forwardSA;
	sdsl::int_vector<> *forwardSA;
	double Cache_percentage;
	int *r2f;
	int *f2r;
	bwi_out *s_for;  //data structure used to manage the FM-index of the (forward) reference string
    bwi_out *s_rev;  //data structure used to manage the FM-index of the reverse reference string

    int ilog2n_reference;  //number of bits needed to represent a reference position

    unordered_map<char, unordered_map<int,int>> reverseOccs;  //caches containing the number of occurrences (value of the internal map)
    unordered_map<char, unordered_map<int,int>> forwardOccs;  //of a character (key of the external map) till a position (key of the internal map)
    														  //these caches are needed to speedup the factorization
    unordered_map<int,int> textPositionsCache;  //caches containing the number of occurrences (value of the internal map)
    //unordered_map<int64_t,Factorization*> factorizations;
    vector<Factorization*> factorizations;
    int maximumFactorsNumber;  //maximum factors number in a factorization
    //vector<Block*> blocks;

    Header header;
    //int *factorizationStarts;  //array containing the single factorization starts positions within the compressed bytestream
    BPlusTree *reverseTree;  //indexes factors by suffixArrayPosition of their reverse in the reverse index, supporting range queries
    BPlusTree *forwardTree; //indexes factors by their suffixArrayPosition in the straight index
    BPlusTree *posTree; //indexes factors by their positions in the reference text

    boost::dynamic_bitset<> charMap; //bitmap used to identify the characters occurring in T
    int lmax;  //used during the factorization

    int cRemappings[256];  //-1 if the character doesn't appear within the text,
        				   //otherwise its code in the compact alphabet of the extra-characters
    vector<int> cInverseRemappings;
    int compactAlphabetSize; //number of characters appearing in the text

    //int textLength; //length of the factorized text
    int bs; //block size

    EncryptedFileBitReader bitReader;
    Portfolio *userPortfolio;


    //The following mutexes are used during the multi-threading factorization process
    mutex charMapLmaxMutex;
    mutex reverseTreeMutex;
    mutex forwardTreeMutex;
    mutex posTreeMutex;
    mutex factorizationsMutex;
    mutex nextToFactorizeMutex;
    mutex coutMutex;

    //Used during the multi-threading factorization process
    int numberOfIndividuals;
    Individual **individuals;
    string *fastaFileNames;
    int nextToFactorize;



};

} /* namespace std */

#endif /* LZINDEX_H_ */
