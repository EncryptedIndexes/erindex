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
#include <map>
#include "Individual.h"
#include <openssl/rsa.h>
#include <openssl/pem.h>

namespace std {

class Portfolio;


class Key{
public:
	uint8_t *getEncryptedValue();
	uint8_t *getClearValue();
	Key(Portfolio *portfolio, uint8_t *encryptedValue);
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
	virtual ~IndividualKey();
	Individual *getIndividual();
	friend class Portfolio;
private:
	Individual *individual;
};


class Portfolio {
public:
	Portfolio(uint64_t userId,RSA *privateKey);
	virtual ~Portfolio();
	friend class Database;
	friend class Key;
	Key *getSystemKey();
	IndividualKey *getIndividualKey(int64_t individualId);
private:
	uint64_t userId;
	RSA *privateKey;
	Key* systemKey;
	map<int64_t,IndividualKey*> individualKeys;  //vector containing the individuals accessible from the user
};






} /* namespace std */

#endif /* PORTFOLIO_H_ */
