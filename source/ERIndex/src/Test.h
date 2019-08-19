/*
 * Test.h
 *
 *  Created on: 05/nov/2014
 *      Author: fernando
 */

#ifndef TEST_H_
#define TEST_H_

#include <vector>
#include <string>

namespace std {

typedef struct PatternSet{
	int length;
	vector<string> patterns;
} PatternSet;

typedef struct TestGroup{
	string inputFileType;
	string databaseRoot;
	string sequencesDirectory;
	string referenceIdentifier;
	string indexFileName;
	bool buildIndex;
	bool naiveSearch;
	bool loadReferenceInMemory;
	bool loadIndexInMemory;
	int blockSize;
	int64_t userId;
	string userPrivateKey;
	string userPassword;
	vector<PatternSet> patternSets;
} TestGroup;

class Test {
public:
	Test(string xmlFilePath);
	virtual ~Test();
	vector<TestGroup> testGroups;
private:
	void loadParameters(string xmlFilePath);
};

} /* namespace std */

#endif /* TEST_H_ */
