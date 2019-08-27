/*
 * Database.h
 *
 *  Created on: 21/nov/2014
 *      Author: fernando
 */

#ifndef DATABASE_H_
#define DATABASE_H_

#include <unordered_map>
#include <map>
#include <vector>
#include "Individual.h"
#include "ReferenceSequence.h"
#include "Portfolio.h"
#include "Catalog.h"

namespace std {

class LZIndex;

class Database {
public:
	Database(string rootDirectory);
	virtual ~Database();
	void login(int64_t userId,string privateKeyFileName);
	Portfolio *getPortfolio();




	void buildIndex(string indexFileName,string sequencesDirectory,int referenceId,int blockSize);
	LZIndex *openIndex(string indexFileName);
	bool verifyIndex(string indexFileName,string sequencesDirectory,int referenceId);

	ReferenceSequence *addReferenceSequence(int id,string referenceSequenceFastaFilePath);
	ReferenceSequence *getReferenceSequence(int id);
	Individual *addIndividual(string &code,string *individualKeyFilePath);
	Individual *getIndividual(int id);
	friend class LZIndex;
private:
	//void loadCatalog();
	Portfolio *getPortfolio(int64_t userId);
	inline bool fileExists (const string& filePath,uint8_t **buffer,bool load);
	void saveCorrespondenceArrayToDisk(string fileName,int *a,int n);
	void buildSuffixArrayCorrespondenceFiles(int n,string reverseSaFileName,string saFileName,string r2fFileName,string f2rFileName);

	string rootDirectory;
	Catalog* catalog;
	uint64_t* loggedUserId;
	//The following two fields were moved to the Catalog class
		//map<int64_t,Individual*> individuals;
		//unordered_map<int64_t,ReferenceSequence> references;




};

} /* namespace std */

#endif /* DATABASE_H_ */
