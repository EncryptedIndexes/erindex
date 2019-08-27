/*
 * Catalog.h
 *
 *  Created on: Aug 19, 2019
 *      Author: fernando
 */

#ifndef CATALOG_H_
#define CATALOG_H_
#include <string>
#include <unordered_map>
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem.hpp>
#include "Individual.h"
#include "ReferenceSequence.h"


namespace std {

class Catalog {
public:
	Catalog(string rootDirectory);
	virtual ~Catalog();
	void addIndividual(Individual *individual);
	void addReferenceSequence(ReferenceSequence *referenceSequence);
	void save();

	map<int64_t,Individual*> individuals;
	unordered_map<int64_t,ReferenceSequence*> references;
private:
	void load();
	void create();
	string rootDirectory;
	string catalogFilePath;
	boost::property_tree::ptree catalogTree;
};

}
#endif /* CATALOG_H_ */