/*
 * ReferenceSequence.h
 *
 *  Created on: 22/nov/2014
 *      Author: fernando
 */

#ifndef REFERENCESEQUENCE_H_
#define REFERENCESEQUENCE_H_
#include <string>
namespace std {

class ReferenceSequence {
public:
	ReferenceSequence();
	virtual ~ReferenceSequence();
	friend class Database;
	friend class LZIndex;
private:
	int id;
	string forwardIndexFileName;
	string reverseIndexFileName;
	string correspondenceFilesName;
	string forwardSuffixArrayFileName;  // usato, solo per il momento, per rimediare alla scarsa efficienza di FM-INDEX v1 nel marcare le righe
};

} /* namespace std */

#endif /* REFERENCESEQUENCE_H_ */
