/*
 * EncryptedFileBitWriter.h
 *
 *  Created on: 27/nov/2014
 *      Author: fernando
 */

#ifndef ENCRYPTEDFILEBITWRITER_H_
#define ENCRYPTEDFILEBITWRITER_H_

#include "FileBitWriter.h"
#include "EncryptionContext.h"
#include "../srcSalsa20/ecrypt-portable.h"
#include "../srcSalsa20/ecrypt-sync.h"



namespace std {


class EncryptedFileBitWriter: public FileBitWriter {
public:
	/* Needed to have a ofstream class member */
	EncryptedFileBitWriter(const EncryptedFileBitWriter&) = delete;
    EncryptedFileBitWriter& operator=(const EncryptedFileBitWriter&) = delete;

    EncryptedFileBitWriter(string fileName, EncryptionContext *ec);
    EncryptionContext *getEncryptionContext();
    void setEncryptionContext(EncryptionContext *ec);

    virtual void putByteIntoBlockBuffer(int ch);
	virtual ~EncryptedFileBitWriter();

private:
	EncryptionContext* eContext;
	void ECRYPT_ctrsetup(ECRYPT_ctx *x,const u8 *ctr);
	void intToBigEndian(uint32_t n, uint8_t *bs);
	void populateKeyStreamBuffer();
};

} /* namespace std */

#endif /* ENCRYPTEDFILEBITWRITER_H_ */
