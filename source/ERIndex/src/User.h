/*
 * User.h
 *
 *  Created on: Aug 28, 2019
 *      Author: fernando
 */

#ifndef USER_H_
#define USER_H_

#include <stdint.h>
#include <string>

namespace std {

class User {
public:
	User();
	virtual ~User();
	friend class Portfolio;
	friend class Database;
	friend class LZIndex;
	uint64_t id; //user internal code
	string username; //username
	string loginValue;  //the string "ERINDEX_LOGIN_VALUE", encrypted with the user's public key, so that
	                    //it can be decrypted only with the user's private key to enforce the login logic
};

} /* namespace std */

#endif /* USER_H_ */
