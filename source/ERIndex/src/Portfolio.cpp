/*
 * Portfolio.cpp
 *
 *  Created on: 21/nov/2014
 *      Author: fernando
 */

#include "Portfolio.h"
#include <openssl/rsa.h>
#include <iostream>


namespace std {

Portfolio::Portfolio(uint64_t userId,RSA *privateKey) {
	this->userId=userId;
	this->privateKey=privateKey;
}

Portfolio::~Portfolio() {
	// TODO Auto-generated destructor stub
}


Key *Portfolio::getSystemKey(){
	return systemKey;
}


IndividualKey *Portfolio::getIndividualKey(int64_t individualId){
	if (individualKeys.count(individualId)>0)
		return individualKeys[individualId];
	else
		return NULL;
}


IndividualKey::IndividualKey(Portfolio *portfolio,Individual *individual,uint8_t *encryptedValue):Key(portfolio,encryptedValue){
	this->individual=individual;
}

IndividualKey::~IndividualKey() {
	// TODO Auto-generated destructor stub
}

Key::Key(Portfolio *portfolio,uint8_t *encryptedValue){
	this->portfolio=portfolio;
	this->encryptedValue=encryptedValue;
	clearValue=NULL;
}

Key::~Key() {
	// TODO Auto-generated destructor stub
}

uint8_t *Key::getEncryptedValue(){
	return encryptedValue;
}

uint8_t *Key::getClearValue(){
	if (clearValue==NULL){
		char *decryptBuffer = new char[64];
		int decryptedKeyLength= RSA_private_decrypt(256, encryptedValue,
								 (unsigned char*)decryptBuffer,portfolio->privateKey , RSA_PKCS1_PADDING);
		clearValue=(uint8_t*)decryptBuffer;
	}
	return clearValue;
}


} /* namespace std */
