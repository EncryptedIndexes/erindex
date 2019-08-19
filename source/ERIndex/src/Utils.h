/*
 * Utils.h
 *
 *  Created on: 24/nov/2014
 *      Author: fernando
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <string>
namespace std {

class Utils {
public:
	Utils();
	virtual ~Utils();
	static string loadFasta(string fastaFileName);
};

} /* namespace std */

#endif /* UTILS_H_ */
