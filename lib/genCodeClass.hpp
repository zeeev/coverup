//
//  genCodeClass.hpp
//  gencodeLib
//
//  Created by Zev Kronenberg on 10/1/15.
//  Copyright © 2015 Zev Kronenberg. All rights reserved.
//

#ifndef genCodeClass_hpp
#define genCodeClass_hpp

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "split.hpp"

class gene;
class transcript;
class exon;
class cds;

class gene{
friend class gcClass;
protected:
  bool proteinCoding;
  int  start;
  int  end;
  int  length;
  int  childrenLen;
  char strand;
  std::string type;
  std::string chr;
  std::string id;
  std::string geneName;
  std::map<std::string, gene *> children;
  std::map<std::string, std::string> attributes;
public:
  virtual bool parseFeature(std::string &);
  virtual void printFeature(void);
  bool isProteinCoding(void);
  bool getLongestChild(gene **);
  bool getChildren(std::vector<gene *> &);
  int getChildrenLength(void);
  int getStart(void);
  int getEnd(void);
  int getLength(void);
  std::string getSeqid(void);
  std::string getGeneName(void);
  
  ~gene(void);
  gene(void);
  
};

class transcript : public gene {
protected:
public:
    bool parseFeature(std::string &);
};


class exon : public transcript {
private:
public:
    bool parseFeature(std::string &);
};

class cds : public exon {
    exon * parent;
};

class indexEntry{
    friend class gcClass;
    std::string id;
    int start     ;
    int end       ;
};

class gcClass{
private:
    bool indexLoaded;
    std::string filename;
    std::string indexname;
    std::ifstream myfile;
    std::vector<indexEntry *> indexDat;
    std::vector<indexEntry *>::iterator ii;
public:
    gcClass();
    gcClass(std::string);
    ~gcClass();
    void index(void);
    bool checkForIndex(void);
    bool loadIndex(void);
    bool getNextGene(gene &);
};


#endif /* genCodeClass_hpp */
