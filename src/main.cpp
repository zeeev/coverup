//
//  main.cpp
//  gencodeLib
//
//  Created by Zev Kronenberg on 10/1/15.
//  Copyright © 2015 Zev Kronenberg. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "genCodeClass.hpp"


int main(int argc, const char * argv[]) {
 

    
    gcClass gtf ("data/gencode.v19.annotation.gtf");

    std::cout << "indexing file\n";    
    gtf.index();
    std::cout << "file was indexed\n";

    gtf.loadIndex();

    
    while(true){
        gene gf;
        gene * transcript;
        if(gtf.getNextGene(gf)){
            if(!gf.isProteinCoding()){
                continue;
            }
            
        gf.printFeature();
        gf.getLongestChild(&transcript);
    
            transcript->printFeature();
    
        }
        else{
            break;
        }
    }

        
    // insert code here...
    
    return 0;
}
