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

	ReferenceSequence *getReferenceSequence(int id);
	Individual *getIndividual(int id);
	friend class LZIndex;
private:
	void loadCatalog();
	Portfolio *getPortfolio(int64_t userId);
	inline bool fileExists (const string& filePath,uint8_t **buffer,bool load);
	string rootDirectory;
	uint64_t* loggedUserId;
	map<int64_t,Individual*> individuals;
	unordered_map<int64_t,ReferenceSequence> references;
};

} /* namespace std */

#endif /* DATABASE_H_ */
