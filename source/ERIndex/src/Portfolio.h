/*
 * Portfolio.h
 *
 *  Created on: 21/nov/2014
 *      Author: fernando
 */

#ifndef PORTFOLIO_H_
#define PORTFOLIO_H_
#include <vector>
#include <string>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include "Individual.h"
#include <openssl/rsa.h>
#include <openssl/pem.h>

namespace std {

class Portfolio;


class Key{
public:
	void setEncryptedValue(uint8_t *encryptedValue);
	uint8_t *getEncryptedValue();
	uint8_t *getClearValue();
	Key(Portfolio *portfolio, uint8_t *encryptedValue);
	Key(Portfolio *portfolio);
	virtual ~Key();
	friend class Portfolio;
protected:
	Portfolio *portfolio;
	uint8_t *encryptedValue; //the individual encryption key, cyphered with the actual user public key
	uint8_t *clearValue; //the individual encryption key in clear

};

class IndividualKey: public Key {
public:
	IndividualKey(Portfolio *portfolio,Individual *individual, uint8_t *encryptedValue);
	IndividualKey(Portfolio *portfolio,Individual *individual);
	virtual ~IndividualKey();
	Individual *getIndividual();
	friend class Portfolio;
private:
	Individual *individual;
};


class Portfolio {
public:
	Portfolio(uint64_t userId,RSA *privateKey);
	Portfolio(uint64_t userId);
	virtual ~Portfolio();
	friend class Database;
	friend class Key;
	void setPublicKey(RSA *publicKey);
	Key *getSystemKey();
	void addSystemKey(uint8_t *clearValue);
	void save(string directory);
	IndividualKey *getIndividualKey(int64_t individualId);
	void addIndividualKey(Individual *individual,uint8_t *clearValue);
	void saveIndividualKey(string directory,uint64_t individualId);
private:
	uint64_t userId;
	RSA *privateKey;
	RSA *publicKey;
	Key* systemKey;
	map<int64_t,IndividualKey*> individualKeys;  //vector containing the individuals accessible from the user
	void saveSystemKey(string directory);
	void saveIndividualKeys(string directory);
};






} /* namespace std */

#endif /* PORTFOLIO_H_ */
