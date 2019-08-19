/*
 * EncryptionContext.h
 *
 *  Created on: 27/nov/2014
 *      Author: fernando
 */

#ifndef ENCRYPTIONCONTEXT_H_
#define ENCRYPTIONCONTEXT_H_
#include "../srcSalsa20/ecrypt-portable.h"
#include <stdint.h>
namespace std {

struct EncryptionContext {
public:
	EncryptionContext(u8 *encKey,uint64_t nonce);
	EncryptionContext(const EncryptionContext &obj);
	virtual ~EncryptionContext();
	u8 *key; //encryption key
	uint64_t nonce; //nonce (used as initialization vector)
	uint32_t ksBlockCounter; //current block index (used to initialize the Salsa20/20 ctr parameter)
	u8 *ksBlockBuffer;  //keystream block buffer: contains the 64 bytes generated by the last Salsa20/20 invocation
	int8_t ksBlockOffset; //current offset of the keystream block buffer
};

} /* namespace std */

#endif /* ENCRYPTIONCONTEXT_H_ */