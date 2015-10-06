//
//  main.cpp
//  
//
//  Created by Zev Kronenberg on 10/1/15.
//  Copyright © 2015 Zev Kronenberg. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <bitset>
#include <vector>
#include "split.hpp"
#include "yagbv.hpp"

#include "genCodeClass.hpp"
#include "tabixpp/tabix.hpp"

struct opts{
  int individuals;
}globalOpts;

void destoryGenomicRange(std::vector<gr *> & grt){
  for(std::vector<gr *>::iterator it = grt.begin(); it != grt.end(); it++){
    delete (*it);
  }
}

void turnOnRanges(std::vector<gr *> & iTranscript, 
		  std::vector<gr *> & iExon,
		  gene * transcript                ){
  
  std::vector<gene *> exons;
  transcript->getChildren(exons);
  
  for(int i = 0 ; i < globalOpts.individuals; i++){
    gr * t = new gr(transcript->getStart(), transcript->getEnd());
    gr * e = new gr(transcript->getStart(), transcript->getEnd());
    t->setRange(0, transcript->getLength());
    iTranscript.push_back(t)                                     ;
    iExon.push_back(e)                                           ;
  } 
  for(std::vector<gene *>::iterator it = exons.begin(); it != exons.end(); it++){
    for(int i = 0 ; i < globalOpts.individuals; i++){
      iExon[i]->setRange((*it)->getStart() - transcript->getStart(), (*it)->getLength());
    }
  }
}

void turnOffRanges(Tabix & tbx, 
		   std::vector<gr *> iTranscript,
		   std::vector<gr *> iExon,
		   gene * transcript ){
  
  stringstream regionss;
  regionss << 1 << ":" << transcript->getStart() << "-" << transcript->getEnd() ;
  string region = regionss.str();
  
  tbx.setRegion(region);
  
  std::string line;
  while(tbx.getNextLine(line)){
    
    std::vector<std::string> linDat = split(line, "\t");
    
    std::map<std::string, std::string> attributes;

    std::vector<std::string> tmp = split(linDat[7], ";");
    for(std::vector<std::string>::iterator iz = tmp.begin(); iz != tmp.end(); iz++){
      std::vector<std::string> kv = split(*iz, "=");
      attributes[ kv[0] ] = kv[1];
    }
    
    std::cerr << "parsed VCF line" << std::endl;
    
    if(attributes.find("END") == attributes.end()){
      continue;
    }
    
    int svstart = atoi(linDat[1].c_str());
    int svend = atoi(attributes["END"].c_str());
    
    // feature contains the sv
    if(svstart >= transcript->getStart() &&
	    svend <= transcript->getEnd()){
std::cerr << "feature contains the sv" << std::endl;
      int offset = svstart - transcript->getStart();
      int svlen  = svend - svstart;

//      transcript->printFeature();
//      std::cerr << "sv start: " << svstart << std::endl;
//      std::cerr << "sv end  : " << svend   << std::endl;
//
//      std::cerr << "offset: " << offset << std::endl;
//      std::cerr << "svlen: " << offset << std::endl;

      for(int i= 9; i < globalOpts.individuals+9; i++){
	if(linDat[i].compare("0|0") != 0){
	  iTranscript[i]->clearRange(offset, svlen);
	  iExon[i]->clearRange(offset, svlen);
	}
      }
      std::cerr << "feature contains the sv OK" << std::endl;    
    }
    // sv contains feature
    else if(svstart <= transcript->getStart() 
	    && svend >= transcript->getEnd()){
       std::cerr << "sv contains feature" << std::endl;

      int len = svend - transcript->getStart();

      transcript->printFeature();
      std::cerr << "sv start: " << svstart << std::endl;
      std::cerr << "sv end  : " << svend   << std::endl;


      for(int i= 9; i < globalOpts.individuals+9; i++){
	if(linDat[i].compare("0|0") != 0){
	  iTranscript[i]->clearRange(0, transcript->getLength());
	  iExon[i]->clearRange(0, transcript->getLength());
	}
      }
    }
    // sv overlaps feature on left
    else if(svstart <= transcript->getStart() 
	    && svend <= transcript->getEnd()){
       std::cerr << "sv overlaps feature on left" << std::endl;
      int len = svend - transcript->getStart();
      for(int i= 9; i < globalOpts.individuals+9; i++){
	if(linDat[i].compare("0|0") != 0){
	  iTranscript[i]->clearRange(0, len);
	  iExon[i]->clearRange(0, len);
	}
      }
    }

    // sv overlaps feature on the right
    
    else{
      int offset = svstart - transcript->getStart();
       std::cerr << "sv overlaps feature on right" << std::endl;
      for(int i= 0; i < globalOpts.individuals; i++){
	if(linDat[i].compare("0|0") != 0){
	  iTranscript[i]->clearRange(offset, transcript->getLength());
	  iExon[i]->clearRange(offset, transcript->getLength());
	}
      }
    }
  }
    std::cerr << "returning from parsed VCF line" << std::endl;
}

int main(int argc, const char * argv[]) {

  globalOpts.individuals = 2504;

  std::string tbf = "data/ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz";
 
  Tabix tbx(tbf);
  
  gcClass gtf ("data/gencode.v23.annotation.gtf");
  
  gtf.index();
  gtf.loadIndex();
    
  bool flag = true;
  
  while(true){
    gene gf;
    gene * transcript;
    flag = gtf.getNextGene(gf);
    if(!flag){
      break;
    }
    if(!gf.isProteinCoding()){
      continue;
    }
    gf.getLongestChild(&transcript);
      
    if(gf.getSeqid().compare("chr1") != 0){
      break;
    }

    std::vector<gr *> iTranscript;
    std::vector<gr *> iExon      ;
    
    turnOnRanges(iTranscript, iExon,       transcript);
    turnOffRanges(tbx, iTranscript, iExon, transcript);
    
    for(int i = 0; i < globalOpts.individuals; i++){
      std::cerr << iTranscript[i]->countOff() << " " << iExon[i]->countOff() << std::endl;
    }
    destoryGenomicRange(iTranscript);
    destoryGenomicRange(iExon);   
  }

  std::cerr << "INFO: coverup finished normally! " << std::endl;
  return 0;
}
