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
#include "User.h"

namespace std {

class LZIndex;

class Database {
public:
	Database(string rootDirectory);
	virtual ~Database();
	static void initialize(string rootDirectory,string systemKeyFile);
	void login(int64_t userId,string privateKeyFileName);
	Portfolio *getPortfolio();

	void buildIndex(string indexFileName,string sequencesDirectory,int referenceId,int blockSize);
	LZIndex *openIndex(string indexFileName);
	bool verifyIndex(string indexFileName,string sequencesDirectory,int referenceId);

	ReferenceSequence *addReferenceSequence(int id,string referenceSequenceFastaFilePath);
	ReferenceSequence *getReferenceSequence(int id);
	Individual *addIndividual(string &code,string *individualKeyFilePath);
	Individual *getIndividual(int id);
	User *addUser(string username,int *userId,string *privateKeyFilePath,string *publicKeyFilePath);
	User *getUser(int userId);
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
	User *loggedUser;


};

} /* namespace std */

#endif /* DATABASE_H_ */
