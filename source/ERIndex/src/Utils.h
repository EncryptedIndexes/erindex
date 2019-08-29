/*
 * Utils.h
 *
 *  Created on: 24/nov/2014
 *      Author: fernando
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <string>
#include "../srcSalsa20/ecrypt-portable.h"
extern "C" {
	#include "../srcFMIndex/common.h"
}



namespace std {


class Utils {
public:
	Utils();
	virtual ~Utils();
	static int fbit_read_simple(FILE *f,int* Bit_buffer_size, uint32* Bit_buffer,int n);
	static int fbit_read24_simple(FILE *f, int* Bit_buffer_size, uint32* Bit_buffer, int n);
	static char *toCString(string s);
	static string loadFasta(string fastaFileName);
	static void generateSalsa20KeyFile(string &keyFilePath);
	static int generateRSAKeyPair(string &privateKeyFileName,string &publicKeyFileName);
};

} /* namespace std */

#endif /* UTILS_H_ */
