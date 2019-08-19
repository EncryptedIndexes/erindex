/*
 * Database.cpp
 *
 *  Created on: 21/nov/2014
 */

#include "Database.h"
#include "Utils.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include "LZIndex.h"
namespace std {


Database::Database(string rootDirectory) {
	this->rootDirectory=rootDirectory;
	this->loggedUserId=NULL;
	loadCatalog();
}

Database::~Database() {
	// TODO Auto-generated destructor stub
}

void Database::login(int64_t userId,string privateKeyFileName){
	//TODO: login logic

	//set the current user
	loggedUserId=new uint64_t;
	*loggedUserId=userId;
}

/*
void Database::buildIndex(string indexFileName,string sequencesDirectory,int referenceId,int blockSize){
	LZIndex index(*this,referenceId,blockSize);
	for( auto it = individuals.begin(), end = individuals.end(); it != end; ++it){
	    Individual *individual=it->second;
	    cout << "Adding individual "<< individual->id << "(" << individual->code << "): ";
	    string sequenceFilePath=sequencesDirectory + "/" + individual->code + "_" + to_string(referenceId) +".fa";
	    Factorization *f=index.addToIndex(individual,sequenceFilePath,1);
	    cout << "\t" << f->factors.size() << " factors" << endl;
	}
	index.save(rootDirectory + "/indexes/"+ indexFileName);
	index.close();
}
*/

void Database::buildIndex(string indexFileName,string sequencesDirectory,int referenceId,int blockSize){
	LZIndex index(*this,referenceId,blockSize);
	Individual **individs=new Individual*[individuals.size()];
	string *fastaFileNames=new string[individuals.size()];
	int i=0;
	for( auto it = individuals.begin(), end = individuals.end(); it != end; ++it){
	    Individual *individual=it->second;
	    individs[i]=individual;
	    string sequenceFilePath=sequencesDirectory + "/" + individual->code + "_" + to_string(referenceId) +".fa";
	    fastaFileNames[i]=sequenceFilePath;
	    i++;
	}
	index.buildIndex(individuals.size(),individs,fastaFileNames);
	index.save(rootDirectory + "/indexes/"+ indexFileName);
	index.close();
}

LZIndex *Database::openIndex(string indexFileName){
	LZIndex *index=new LZIndex(*this);
	index->open(rootDirectory + "/indexes/"+ indexFileName);

	return index;
}


bool Database::verifyIndex(string indexFileName,string sequencesDirectory,int referenceId){
	LZIndex *index=new LZIndex(*this);
	index->open(rootDirectory + "/indexes/"+ indexFileName);
	bool ok=true;
	vector<int64_t> ids=index->getIndividualIds();
	for (int i=0;i<ids.size();i++){

		string s=index->getSequence(ids[i]);
		Individual *individual=individuals[ids[i]];
		cout << "Verifying the individual " << i << "(" << individual->code << ")" << endl;
		string sequenceFilePath=sequencesDirectory + "/" + individual->code + "_" + to_string(referenceId) +".fa";
		string t=Utils::loadFasta(sequenceFilePath);
		if (s!=t){
			ok=false;
			int minSize;
			if (s<t)
				minSize=s.size();
			else
				minSize=t.size();
			for (int i=0;i<minSize;i++){
				if (s[i]!=t[i]){
					cout << "ERROR at " << i << endl;
					break;
				}
			}
		}
	}
	return ok;
}


void Database::loadCatalog(){
	string catalogFilePath=rootDirectory+"/catalog.xml";
	boost::property_tree::ptree pTree;
	boost::property_tree::read_xml(catalogFilePath,pTree);
	//read the test parameters from the XML file
	//Using boost::property_tree
	for( auto const& individualChildNode: pTree.get_child("catalog.individuals") ) {
			boost::property_tree::ptree individualSubTree = individualChildNode.second;
			int id=individualSubTree.get<int>("id");
			string code=individualSubTree.get<string>("code");
			Individual *individual=new Individual(id,code);
			individuals[id]=individual;
	}
	for( auto const& referenceChildNode: pTree.get_child("catalog.references") ) {
			boost::property_tree::ptree referenceSubTree = referenceChildNode.second;
			ReferenceSequence seq;
			string referencesPath=rootDirectory+ "/references/";
			seq.id=referenceSubTree.get<int>("id");
			seq.forwardIndexFileName=referencesPath+referenceSubTree.get<string>("forwardIndexFileName");
			seq.reverseIndexFileName=referencesPath+referenceSubTree.get<string>("reverseIndexFileName");
			seq.correspondenceFilesName=referencesPath+referenceSubTree.get<string>("correspondenceFilesName");
			seq.forwardSuffixArrayFileName=referencesPath+referenceSubTree.get<string>("forwardSuffixArrayFileName");
			references[seq.id]=seq;
	}
}


Portfolio *Database::getPortfolio (){
	if (loggedUserId!=NULL)
		return getPortfolio(*loggedUserId);
	else
		return NULL;
}


Portfolio *Database::getPortfolio (int64_t userId){
/*
	char * decrypt = (char*)malloc(8);
	int temp = RSA_private_decrypt(256, (unsigned char*)encrypt, (unsigned char*)decrypt,rsa_prikey , RSA_PKCS1_OAEP_PADDING);
	if (!RAND_bytes(key, sizeof key)) {
	    cout << "Aiutooo";
	}
*/
	RSA *privateKey = NULL;
	string filePath=rootDirectory+"/security/user_"+ boost::lexical_cast<std::string>(userId)+ ".pvt";
	FILE *privateKeyFile = fopen(filePath.c_str(),"rb");
	//load user's private key
	if (PEM_read_RSAPrivateKey(privateKeyFile, &privateKey, NULL, NULL) != NULL){
	//if (fileExists(filePath,buffer)){
		Portfolio *portfolio=new Portfolio(userId,privateKey);
		fclose(privateKeyFile);
		//load system key
		filePath=rootDirectory+"/security/portfolio_"+ to_string(userId)+"_s.key";
		uint8_t *buffer=NULL;
		if (fileExists(filePath,&buffer,true)){
			Key *sk=new Key(portfolio,buffer);
			portfolio->systemKey=sk;
		}
		//load individual keys
		for (auto it=individuals.begin(); it!=individuals.end(); ++it){
			//filePath=rootDirectory+"/security/portfolio_"+ boost::lexical_cast<std::string>(userId)+"_"+boost::lexical_cast<std::string>(individuals[i]->id)+ ".key";
			filePath=rootDirectory+"/security/portfolio_"+ to_string(userId)+"_"+to_string(it->first)+ ".key";
			uint8_t *buffer=NULL;
			if (fileExists(filePath,&buffer,true)){
				IndividualKey *ik=new IndividualKey(portfolio,it->second,buffer);
				portfolio->individualKeys[it->second->id]=ik;
			}
		}
		return portfolio;
	} else{
	  fclose(privateKeyFile);
	  throw;
	}
}


inline bool Database::fileExists (const string& filePath,uint8_t **buffer,bool load) {
	ifstream f(filePath.c_str());
	if (f.good()) {
		if (load){
			f.seekg (0, f.end);
			int fileSize = f.tellg();
			*buffer=new uint8_t[fileSize];
			f.seekg (0, f.beg);
			f.read((char *)*buffer,fileSize);
		}
        f.close();
        return true;
    } else {
        return false;
    }
}


ReferenceSequence *Database::getReferenceSequence(int id){
	if (references.count(id)>0)
		return &references[id];
	else
		return NULL;
}


Individual *Database::getIndividual(int id){
	if (individuals.count(id)>0)
		return individuals[id];
	else
		return NULL;
}




} /* namespace std */
