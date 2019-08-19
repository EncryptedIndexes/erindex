/*
 * EncryptedFileBitReader.h
 *
 *  Created on: 27/nov/2014
 *      Author: fernando
 */

#ifndef ENCRYPTEDFILEBITREADER_H_
#define ENCRYPTEDFILEBITREADER_H_

#include "FileBitReader.h"
#include "EncryptionContext.h"
#include "../srcSalsa20/ecrypt-portable.h"
#include "../srcSalsa20/ecrypt-sync.h"
namespace std {

class EncryptedFileBitReader: public FileBitReader {
public:
	virtual ~EncryptedFileBitReader();
	EncryptedFileBitReader();
	EncryptedFileBitReader(string fileName, EncryptionContext *ec);
	EncryptionContext *getEncryptionContext();
	void open(string fileName,EncryptionContext *ec);
	void setEncryptionContext(EncryptionContext *ec);
	virtual int getByteFromBlockBuffer();
private:
	EncryptionContext* eContext;
	void ECRYPT_ctrsetup(ECRYPT_ctx *x,const u8 *ctr);
	void intToBigEndian(uint32_t n, uint8_t *bs);
	void populateKeyStreamBuffer();
};

} /* namespace std */

#endif /* ENCRYPTEDFILEBITREADER_H_ */
