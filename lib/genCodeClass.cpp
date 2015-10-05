//
//  genCodeClass.cpp
//  gencodeLib
//
//  Created by Zev Kronenberg on 10/1/15.
//  Copyright © 2015 Zev Kronenberg. All rights reserved.
//

#include "genCodeClass.hpp"


gene::gene(void){
    proteinCoding = false;
    childrenLen = 0;
}

bool gene::isProteinCoding(void){
    
    if(proteinCoding){
        return true;
    }
    else{
        return false;
    }
    
}


int gene::getChildrenLength(void){
    if(childrenLen > 0){
        return childrenLen;
    }
    else{
        for(std::map<std::string, gene *>::iterator it = children.begin();
            it != children.end(); it++){
            childrenLen += it->second->length;
        }
    }
    return childrenLen;
}

bool gene::getChildren(std::vector<gene *> & v){
    if(children.empty()){
        return false;
    }
    
    for(std::map<std::string, gene *>::iterator it = children.begin();
        it != children.end(); it++){
        v.push_back(it->second);
    }

    return true;
}

bool gene::getLongestChild(gene ** pt){
    if(children.empty()){
        return false;
    }

    int len = -1;
    gene * l = NULL;
    

    for(std::map<std::string, gene *>::iterator it = children.begin();
        it != children.end(); it++){
        if(it->second->length > len){
            l = it->second;
            len = l->length;
        }
    }
    
    *pt = l;
    return true;
}


bool transcript::parseFeature(std::string & line){
    gene::parseFeature(line);
    type = "transcript";
    id   = attributes["transcript_id"];
    return true;
}

bool exon::parseFeature(std::string & line){
    gene::parseFeature(line);
    type = "exon";
    id   = attributes["exon_id"];
    return true;
}

void gene::printFeature(void){
    std::cerr << "type: " << type << std::endl;
    std::cerr << "contig: " << chr << std::endl;
    std::cerr << "start: " << start << std::endl;
    std::cerr << "start: " << end << std::endl;
    std::cerr << "length: " << length << std::endl;
    std::cerr << "strand: " << strand << std::endl;
    std::cerr << "id: "     << id << std::endl;
    std::cerr << "geneName: "     << geneName << std::endl;

    std::cerr << std::endl;
    

    
}


bool gene::parseFeature(std::string & line){
    
    std::vector<std::string> ldat = split(line, '\t');
    
    chr    = ldat[0];
    start  = atoi(ldat[3].c_str());
    end    = atoi(ldat[4].c_str());
    length = end - start;
    type   = "gene";
    strand = ldat[6][0];
    
    std::vector<std::string> firstAttributSplit = split(ldat.back(), ";");
    
    for(std::vector<std::string>::iterator iz = firstAttributSplit.begin();
        iz != firstAttributSplit.end(); iz++){
        
        if((*iz).empty()){
            continue;
        }
        
        std::vector<std::string> kv = split((*iz), ' ');
       
        int firstNonWhite = 0;
        
        while(kv[firstNonWhite].empty()){
            firstNonWhite += 1;
        }
        
        attributes[kv[firstNonWhite]] = kv[(firstNonWhite +1)];
    }
    
    id = attributes["gene_id"];
    geneName = attributes["gene_name"];
    
    return true;
}

gcClass::gcClass(void){
    std::cerr << "FATAL: did not call constructor with GTF" << std::endl;
    exit(1);
}

gene::~gene(void){
    for (std::map<std::string, gene *>::iterator it = children.begin() ;
         it != children.end(); it++) {
        delete  it->second;
    }
}

gcClass::gcClass(std::string file){
    filename = file;
    indexLoaded  = false;
    indexname = filename + ".gindx";
    myfile.open(filename);
}

gcClass::~gcClass(void){
    for(std::vector<indexEntry *>::iterator it = indexDat.begin();
        it != indexDat.end(); it++){
        delete (*it);
    }
    myfile.close();
}

bool gcClass::loadIndex(void){

    if(indexLoaded == true){
        return true;
    }
    
    if(checkForIndex()){

        std::ifstream indexf(indexname);
        std::string line;
        
        if (myfile.is_open())
        {
            while ( getline (indexf,line) )
            {
                indexEntry * datum = new indexEntry;
                std::vector<std::string> ldat = split(line, '\t');
                datum->id    = ldat[0];
                datum->start = atoi(ldat[1].c_str());
                datum->end   = atoi(ldat[2].c_str());
                indexDat.push_back(datum);
                
            }
        }
        std::cerr << "INFO: Loaded GTF index with " << indexDat.size() << " genes." << std::endl;
        indexf.close();
        ii = indexDat.begin();
        
	indexLoaded = true;
	return true;
    }
    else{
        return false;
    }
}

bool gcClass::checkForIndex(void){

    std::ifstream f(indexname);
        if (f.good()) {
            f.close();
            return true;
        } else {
            f.close();
            return false;
        }   
}


void gcClass::index(){
    if(filename.empty()){
        std::cerr << "FATAL: no file to index" << std::endl;
        exit(1);
    }
    
    if(checkForIndex()){
        return;
    }
        
    std::string   line;
    std::ofstream index (indexname);
    
    std::string genename   ;
    long int start  = 0;
    long int end    = 0;
    int lastLineLen = 0;
    
    
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            
            if(line[0] == '#'){
                continue;
            }
            std::vector<std::string> ldat = split(line , '\t');

            if(! (ldat[2].compare("gene") == 0)){
                lastLineLen = myfile.tellg() ;
                lastLineLen -=  (line.size());
                continue;
            }
            std::vector<std::string> info = split(ldat.back(), '"');
         
            long int tmpstart = myfile.tellg() ;
            tmpstart = tmpstart - (line.length()+1);
            
            if(start == 0){
                start    = tmpstart;
                genename = ldat[2];
                continue;
            }
            end = lastLineLen - 1;
            index << genename << "\t" << start << "\t" << end << std::endl;
            genename = ldat[2];
            start = tmpstart;
        }
        index.close();
    }
    else{
        std::cout << "FATAL: Unable to open GTF file" << std::endl;
        exit(1);
    }
}

bool gcClass::getNextGene(gene  & g){
    
    if(ii != indexDat.end()){
        
        std::string line;
        myfile.seekg((*ii)->start, std::ios::beg);

        getline(myfile, line);
        
        g.parseFeature(line);
        
        while(myfile.tellg() != (*ii)->end){
        
            getline(myfile, line);
            
            std::vector<std::string> ldat = split(line , '\t');
            
            gene * unkownF;
            if(ldat[2].compare("transcript") == 0){
                unkownF = new transcript;
                unkownF->parseFeature(line);
                g.children[unkownF->id] = unkownF;
            }
            else if (ldat[2].compare("exon") == 0){
                unkownF = new exon;
                unkownF->parseFeature(line);
                
                // dropping exon in transcript in gene :/
                
                g.children[ unkownF->attributes["transcript_id"] ]->children[unkownF->id] = unkownF;
            }
            else if (ldat[2].compare("CDS") == 0){
                g.proteinCoding = true;
                
                
                unkownF = new cds;
                unkownF->parseFeature(line);
                
                // dropping cds in exon in transcript in gene
            
                g.children[ unkownF->attributes["transcript_id"] ]->children[unkownF->attributes["exon_id"]]->children[unkownF->id] = unkownF;
            }
        }
        
        ii++;
        return true;
    }
    else{
        return false;
    }
}
