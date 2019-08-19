/*
 * Individual.h
 *
 *  Created on: 21/nov/2014
 *      Author: fernando
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <stdint.h>
#include <string>

namespace std {

class Individual {
public:
	Individual();
	Individual(uint64_t id,string code);
	virtual ~Individual();
	friend class Portfolio;
	friend class Database;
	friend class LZIndex;
	uint64_t id; //individual internal code (used to encrypt his genome)
	string code; //mnemonic code of the individual (i.e. "NA06964")
private:

};

} /* namespace std */

#endif /* INDIVIDUAL_H_ */
