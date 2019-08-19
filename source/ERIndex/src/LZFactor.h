/*
 * LZFactor.h
 *
 *  Created on: 21/set/2014
 *      Author: fernando
 */

#ifndef LZFACTOR_H_
#define LZFACTOR_H_

namespace std {

class Factorization;

typedef struct Block;

class LZFactor {
	friend class LZIndex;
	friend Block;
	friend Factorization;
public:
	//used to construct standard factors (substrings of the reference strings)
	LZFactor(Factorization *factorization,int suffixArrayPosition,int length,int textPosition);
	//used to construct extra-factors (characters not existing in the reference string)
	LZFactor(Factorization *factorization,char letter,int textPosition);
	LZFactor(Factorization *factorization);
	string getString();
	int getLength();
	virtual ~LZFactor();

	friend ostream& operator<<(ostream& os, const LZFactor& lz);

private:
	//used for standard factors
	int suffixArrayPosition;
	unsigned int length;
	//used for extra factors
	char letter;
	Factorization *factorization;  //LZIndex* lzIndex;
	//starting position of this factor in the factorized text
	int textPosition;

};

} /* namespace std */



#endif /* LZFACTOR_H_ */
