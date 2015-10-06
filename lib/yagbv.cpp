//
//  yagbv.cpp
//  yagbv
//
//  Created by Zev Kronenberg on 10/5/15.
//  Copyright © 2015 Zev Kronenberg. All rights reserved.
//
//  errors lurk!

#include "yagbv.hpp"


int gr::popcount(uint8_t x)
{
    int c = 0;
    for (; x > 0; x &= x -1) c++;
    return c;
}

gr::~gr(void){
    delete[] datv;
}

gr::gr(int b, int e){
    start = b;
    end   = e;
    
    
    len = end - start;
     
    if(len < 1){
        std::cerr << "FATAL: range must be greater than zero." << std::endl;
        exit(1);
    }
    
    int div  = len / 8;
    int mod  = len % 8;

    offmask = MASK;
    offmask = offmask >> mod;
    
    if(len > 8 && mod > 0){
        div += 1;
    }
    
    if(div == 0){
        div = 1;
    }
    
    alen = div;
    
    datv = new uint8_t [div];
    this->off();
}

void gr::on(void){
    for(int i = 0; i < alen; i++){
        datv[i] = 8;
    }
}

void gr::off(void){
    for(int i = 0; i < alen; i++){
        datv[i] = 0;
    }
}

void gr::set(int i){
   
    datv[i/8] |= 1 << (7-(i%8));
    
}

void gr::clear(int i){
    datv[i/8] &= ~(1 << (7 - i%8));
}

void gr::setRange(int s, int e){
    
    for(int i = s; i <= e; i++){
        set(i);
    }
    
}

void gr::clearRange(int s, int e){
    
    for(int i = s; i <= e; i++){
        clear(i);
    }
    
}


int gr::countOn(void){

    int totalCount = 0;
    
    for(int i = 0; i < alen; i++){
        totalCount += popcount(datv[i]);
    }
    return totalCount;
}


int gr::countOff(void){
    
    int totalCount = 0;
    
    for(int i = 0; i < alen-1; i++){
        totalCount += popcount(~datv[i]);
    }
    
    totalCount += popcount( (~offmask & datv[alen]) );
    
    return totalCount;
}




void gr::print(void){
    for(int i = 0; i < alen; i++){
        std::cerr << int(datv[i]) << " " ;
    }
    std::cerr << std::endl;
}

std::ostream& operator<<(std::ostream& out, gr& g){
    
    for(int i = 0; i < g.alen; i++){
        
        out << int(g.datv[i]);
    }
    return out;
}
