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
#include <boost/filesystem.hpp>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include "LZIndex.h"

extern "C" {
	int createFMIndex(char *infile_name, char *outfile_name,int bs_lev1, int bs_lev2, float mc_freq,
			int save_sa_to_disk);
}

namespace std {


Database::Database(string rootDirectory) {
	this->rootDirectory=rootDirectory;
	this->loggedUserId=NULL;
	this->loggedUser=NULL;

	//Check if the database root directory exists
	if (!boost::filesystem::is_directory(rootDirectory))
		throw runtime_error("The root directory " + rootDirectory + " does not exist");

	//Check if the following sub-directories exist:
	// - references
	// - security
	// - indexes
	// and, if required, create them.
	string subdir=rootDirectory + "/references";
	if (!boost::filesystem::is_directory(subdir))
		if (!boost::filesystem::create_directory(subdir))
			throw runtime_error("Error while creating directory "+subdir);
	subdir=rootDirectory + "/security";
	if (!boost::filesystem::is_directory(subdir))
		if (!boost::filesystem::create_directory(subdir))
			throw runtime_error("Error while creating directory "+subdir);

	string systemKeyPath=rootDirectory+"/security/system.key";
	//Checks if the system key exists and creates it if necessary
	if (!boost::filesystem::is_regular_file(systemKeyPath))
		Utils::generateSalsa20KeyFile(systemKeyPath);

	subdir=rootDirectory + "/indexes";
	if (!boost::filesystem::is_directory(subdir))
		if (!boost::filesystem::create_directory(subdir))
			throw runtime_error("Error while creating directory "+subdir);

	this->catalog= new Catalog(rootDirectory);

	//Test
	//string s="MANU";
	//Individual *ind=addIndividual(s, NULL);
}


void Database::initialize(string rootDirectory,string systemKeyFile){
	//Check if database root directory exists and is empty
	//Ensures that the database root directory exists
	if (!boost::filesystem::is_directory(rootDirectory)){
		cout << "The root directory " + rootDirectory + " does not exist";
		boost::filesystem::create_directories(rootDirectory);
		cout << "-->root directory has been successfully created";
	}
	//Ensures that root directory is empty
	if(!boost::filesystem::is_empty(rootDirectory))
		throw runtime_error("Can't initialize a new database: the root directory " + rootDirectory + " is not empty");

	//Creates the references subdirectory
	string subdir=rootDirectory + "/references";
	if (!boost::filesystem::create_directory(subdir))
		throw runtime_error("Error while creating directory "+subdir);

	//Creates the security subdirectory
	subdir=rootDirectory + "/security";
	if (!boost::filesystem::create_directory(subdir))
		throw runtime_error("Error while creating directory "+subdir);
	//Creates the Salsa20 system key
	string systemKeyDatabasePath=rootDirectory+"/security/system.key";
	if (boost::filesystem::is_regular_file(systemKeyFile))
		boost::filesystem::copy_file(systemKeyFile,systemKeyDatabasePath,boost::filesystem::copy_option::overwrite_if_exists);
	else
		Utils::generateSalsa20KeyFile(systemKeyDatabasePath);

	//Creates the indexes subdirectory
	subdir=rootDirectory + "/indexes";
	if (!boost::filesystem::create_directory(subdir))
		throw runtime_error("Error while creating directory "+subdir);

	//Creates an empty XML database catalog
	Catalog catalog(rootDirectory);

}



Database::~Database() {
	// TODO Auto-generated destructor stub
}


//TODO: replace userId with username
void Database::login(int64_t userId,string privateKeyFileName){

	//set the current user
	loggedUserId=new uint64_t;
	*loggedUserId=userId;

	loggedUser=getUser(userId);
	if (loggedUser==NULL)
		throw new runtime_error("The user "+ to_string(userId)+" does not exist!");

	//TODO: implement login: decrypt the user object login value field and verify that the value
	// decrypted value is exactly ERINDEX_LOGIN_VALUE

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
	Individual **individs=new Individual*[catalog->individuals.size()];
	string *fastaFileNames=new string[catalog->individuals.size()];
	int i=0;
	for( auto it = catalog->individuals.begin(), end = catalog->individuals.end(); it != end; ++it){
	    Individual *individual=it->second;
	    individs[i]=individual;
	    string sequenceFilePath=sequencesDirectory + "/" + individual->code + "_" + to_string(referenceId) +".fa";
	    fastaFileNames[i]=sequenceFilePath;
	    i++;
	}
	index.buildIndex(catalog->individuals.size(),individs,fastaFileNames);
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
		Individual *individual=catalog->individuals[ids[i]];
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


/*
void Database::loadCatalog(){
	//Check if the catalog XML file exists in the root directory.
	//If it doesn't exist, create an empty catalog XML file, otherwise load it
	string catalogFilePath=rootDirectory+"/catalog.xml";

	boost::property_tree::ptree pTree;
	if (boost::filesystem::is_regular_file(catalogFilePath)){

		boost::property_tree::read_xml(catalogFilePath,pTree);

		//read parameters from the XML file
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
	} else {

		boost::property_tree::ptree catalogElement;
		boost::property_tree::ptree referencesElement;
		boost::property_tree::ptree individualsElement;
		catalogElement.push_back(std::make_pair("references", referencesElement));
		catalogElement.push_back(std::make_pair("individuals", individualsElement));

		pTree.push_back(std::make_pair("catalog",catalogElement));

		boost::property_tree::write_xml(catalogFilePath, pTree);
	}
}
*/

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
		for (auto it=catalog->individuals.begin(); it!=catalog->individuals.end(); ++it){
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



void Database::saveCorrespondenceArrayToDisk(string fileName,int *a,int n){
	FileBitWriter bw(fileName);
	bw.open();
	int psize = int_log2(n);
	for(int i=0;i<n;i++)
	    bw.write(psize,a[i]);
	bw.flush();
	bw.close();
}






__inline__ int get_suffix(FILE *Safile, int n,int i)
{
  //int pos,rem;
  long pos,rem;

  //int fbit_read_simple(FILE *,int*, uint32*,int);
  int Pointer_size=int_log2(n);
  uint32* Bit_buffer=(uint32*)malloc(sizeof(uint32));
  int*  Bit_buffer_size=(int*)malloc(sizeof(int));  /* number of unread/unwritten bits in Bit_buffer */

  pos = (long)i*(long)Pointer_size;
  rem = int(pos& 0x7);         // last 3 bits

  fseek(Safile, pos>>3, SEEK_SET);

  *Bit_buffer= (uint32) 0;
  *Bit_buffer_size=0;

  if(rem)
    Utils::fbit_read_simple(Safile,Bit_buffer_size,Bit_buffer,rem);

  int val=Utils::fbit_read_simple(Safile, Bit_buffer_size,Bit_buffer,Pointer_size);
  free(Bit_buffer);
  free(Bit_buffer_size);
  return val;
}




void Database::buildSuffixArrayCorrespondenceFiles(int n,string reverseSaFileName,string saFileName,string r2fFileName,string f2rFileName){

	FILE *revSaFile = fopen(reverseSaFileName.c_str(), "a+b");
	fseek(revSaFile, 0L, SEEK_SET);

	int *saSuffixes=new int[n];
	//fill_n(saSuffixes,n,-1);
	for (int i=0;i<n;i++)
		saSuffixes[i]=-1;

	cout << "Reading suffix array from disk" << endl;
	cout << "It may take a few minutes" << endl;

	FILE *saFile = fopen(saFileName.c_str(), "a+b");
	fseek(saFile, 0L, SEEK_SET);
	for (int i=0;i<n;i++){
		if (i%10000000==0)
			cout << "Reading the "<< (i+1) << "-th suffix array element"<< endl;

		int s=get_suffix(saFile,n,i);
		saSuffixes[s]=i;
	}
	fclose(saFile);

	//Calcolo i due vettori r2f ed f2r, che consentono rispettivamente di passare da un suffisso dellla stringa reverse al corrispondente della straight e viceversa
	int *r2f=new int[n];
	int *f2r=new int[n];
	int revSa_di_i;
	int j;
	cout << "Computing R2F and F2R" << endl;
	for (int i=0;i<n;i++){
		if (i%10000000==0)
			cout << "Computing the "<< (i+1) << "-th R2F and F2R element"<< endl;
		revSa_di_i=get_suffix(revSaFile,n,i);
		if (saSuffixes[n-revSa_di_i-1]==-1){
			cout << "NNCUC" << endl;
			exit(1);
		}
		j=saSuffixes[n-revSa_di_i-1]; //suffix related to the same text position
		r2f[i]=j;
		f2r[j]=i;
	}

	fclose(revSaFile);

	cout << "Saving R2F and F2R to disk" << endl;
	//Salvo i due array su disco
	saveCorrespondenceArrayToDisk(r2fFileName,r2f,n);
	saveCorrespondenceArrayToDisk(f2rFileName,f2r,n);

	free(saSuffixes);
	free(r2f);
	free(f2r);
}




ReferenceSequence *Database::addReferenceSequence(int id,string referenceSequenceFastaFilePath){
	if (catalog->references.count(id)>0)
		throw new runtime_error("Reference sequence with id " +  to_string(id) + "already exists in the database!");

	if (!boost::filesystem::is_regular_file(referenceSequenceFastaFilePath))
		throw new runtime_error("The given reference sequence file is not valid or not accessible!");

	seqan::SeqFileIn seqFileIn(referenceSequenceFastaFilePath.c_str());
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::IupacString> seqs;
	seqan::readRecords(ids, seqs, seqFileIn);
	vector<string> sequenceFileNames;
	vector<string> sequenceDescriptions;
	int size=length(seqs);
	if (size>0){

		//Gets the sequence length
		int n=seqan::length(seqs[0]);

		boost::filesystem::path p(referenceSequenceFastaFilePath);
		string fileName=p.stem().string();
		string indexFilePath=this->rootDirectory + "/references/" + fileName + ".bwi";


		//Creates the FM-index and the suffix array of the sequence and saves them to disk
		//Creates a temporary file containing only the reference sequence characters
		string sequenceDescription(seqan::toCString(ids[0]));
		sequenceDescriptions.push_back(sequenceDescription);
		char *tfn=std::tmpnam(nullptr);


		std::ofstream temporaryFileStream(tfn);
		temporaryFileStream << seqs[0];
		temporaryFileStream.flush();
		temporaryFileStream.close();



		string temporaryFileName=tfn;
		createFMIndex(Utils::toCString(temporaryFileName), Utils::toCString(indexFilePath),
				8192, 1024, 0.02,1);

		//Creates the FM-index and the suffix array of the sequence reverse complement and saves them to disk

		tfn=std::tmpnam(nullptr);
		std::ofstream temporaryFileStream2(tfn);
		seqan::IupacString revSeq=seqs[0];
		reverseComplement(revSeq);
		temporaryFileStream2 << revSeq;
		temporaryFileStream2.flush();
		temporaryFileStream2.close();

		string revIndexFilePath=this->rootDirectory + "/references/" + fileName + "_rev.bwi";
		temporaryFileName=tfn;
		createFMIndex(Utils::toCString(temporaryFileName), Utils::toCString(revIndexFilePath),
						8192, 1024, 0.02,1);

		//builds the R2F and F2R suffix array correspondence files
		string revSaFileName=this->rootDirectory + "/references/" + fileName + "_rev.sa";
		string saFileName=this->rootDirectory + "/references/" + fileName + ".sa";
		string r2fFileName=this->rootDirectory + "/references/" + fileName + ".r2f";
		string f2rFileName=this->rootDirectory + "/references/" + fileName + ".f2r";
		buildSuffixArrayCorrespondenceFiles(n, revSaFileName, saFileName, r2fFileName, f2rFileName);

		//deletes from disk the reverse suffix array file
		boost::filesystem::remove(revSaFileName);

		//adds the reference sequence to the catalog
		ReferenceSequence *ref=new ReferenceSequence();
		ref->id=id;
		ref->forwardIndexFileName=indexFilePath;
		ref->reverseIndexFileName=revIndexFilePath;
		ref->correspondenceFilesName=this->rootDirectory + "/references/" + fileName;
		ref->forwardSuffixArrayFileName=saFileName;
		catalog->addReferenceSequence(ref);
		catalog->save();

	}

}



ReferenceSequence *Database::getReferenceSequence(int id){
	if (catalog->references.count(id)>0)
		return catalog->references[id];
	else
		return NULL;
}


Individual *Database::getIndividual(int id){
	if (catalog->individuals.count(id)>0)
		return catalog->individuals[id];
	else
		return NULL;
}

Individual *Database::addIndividual(string &code, string *individualKeyFilePath){
	//verify if the individual already exists and retrieve the greatest id in the catalog
	int greatestId=0;
	for (std::map<int64_t,Individual*>::iterator it=catalog->individuals.begin(); it!=catalog->individuals.end(); ++it){

		if (it->second->code == code)
			throw new runtime_error("An individual with code " +  code + "already exists in the database!");

		if (it->second->id>greatestId)
			greatestId=it->second->id;
	}
	int newId=greatestId+1;

	//creates, if needed, the Salsa 20 key for the individual
	string *ifp = new string(rootDirectory + "/security/" + code + ".key");
	if (individualKeyFilePath ==NULL){
		Utils::generateSalsa20KeyFile(*ifp);
	}
	else{
		if (!boost::filesystem::is_regular_file(*individualKeyFilePath))
				throw new runtime_error("The given key file path for the individual "+ code +  " is not valid or not accessible!");
		boost::filesystem::copy_file(*individualKeyFilePath,*ifp,boost::filesystem::copy_option::overwrite_if_exists);
	}

	//adds the individual to the catalog
	Individual *ind=new Individual(newId,code);
	catalog->addIndividual(ind);
	catalog->save();

}






User *Database::addUser(string username,int *userId,string *privateKeyFilePath,string *publicKeyFilePath){
	User *u=new User();
	u->username=username;
	//if the parameter userid was given, check if it exists in the database
	if (userId!=NULL){
		if (catalog->users.count(*userId)>0)
			throw new runtime_error("The user with id " +  to_string(*userId) + "already exists in the database!");
		u->id=*userId;
	}
	else{
		int greatestId=0;
		for (std::map<int64_t,User*>::iterator it=catalog->users.begin(); it!=catalog->users.end(); ++it){

			if (it->second->username == username)
				throw new runtime_error("The user " +  username + "already exists in the database!");

			if (it->second->id>greatestId)
				greatestId=it->second->id;
		}
		u->id=greatestId+1;
	}




	string dbPrivateKeyFilePath=this->rootDirectory + "/security/user_"+   to_string(u->id) + ".pvt";
	string dbPublicKeyFilePath=this->rootDirectory + "/security/user_"+ to_string(u->id) + ".pub";
	if (privateKeyFilePath == NULL || publicKeyFilePath == NULL)
		Utils::generateRSAKeyPair(dbPrivateKeyFilePath, dbPublicKeyFilePath);
	else{
		boost::filesystem::copy_file(*privateKeyFilePath,dbPrivateKeyFilePath,boost::filesystem::copy_option::overwrite_if_exists);
		boost::filesystem::copy_file(*publicKeyFilePath,dbPublicKeyFilePath,boost::filesystem::copy_option::overwrite_if_exists);
	}

	//TODO:encrypts the value ERINDEX_LOGIN_VALUE with the user's public key
	//and stores it into the loginValue field

	catalog->addUser(u);
	catalog->save();


}


User *Database::getUser(int userId){
	if (catalog->users.count(userId)>0)
		return catalog->users[userId];
	else
		return NULL;
}





} /* namespace std */

