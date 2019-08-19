/*
 * Utils.cpp
 *
 *  Created on: 24/nov/2014
 *      Author: fernando
 */

#include "Utils.h"
#include <string>
#include <iostream>
#include <fstream>
namespace std {

Utils::Utils() {
	// TODO Auto-generated constructor stub

}

Utils::~Utils() {
	// TODO Auto-generated destructor stub
}

/**
 * Loads in a buffer of size sequenceLength + k-1 the sequence contained in the fastaFileName,
 * so that the buffer contains all the characters to compute the initial k-mers of the last text suffixes
 * @param fastaFileName
 * @throws Exception
 */
string Utils::loadFasta(string fastaFileName) {
	std::ifstream in(fastaFileName, std::ios::in);
	if (in) {
		std::string sequence,line;
		std::getline( in, line );
		while( std::getline( in, line ).good() ){
		        sequence += line;
		}
		return sequence;
	}
}

} /* namespace std */
