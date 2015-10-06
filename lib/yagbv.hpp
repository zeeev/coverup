//
//  yagbv.hpp
//  yagbv
//
//  Created by Zev Kronenberg on 10/5/15.
//  Copyright © 2015 Zev Kronenberg. All rights reserved.
//

#ifndef yagbv_hpp
#define yagbv_hpp

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <ostream>

#define MASK 255

class gr{

private:
    int start       ;
    int end         ;
    int len         ;
    int alen        ;
    uint8_t   offmask ; /* the mask for the trailing */
    uint8_t * datv  ;
    void set(int )  ;
    void clear(int) ;
    int popcount(uint8_t);
public:
    gr(int, int);
    ~gr(void);
    int countOn(void);
    int countOff(void);
    
    void on(void );
    void off(void);
    void print(void);
    void setRange(int, int);
    void clearRange(int, int);
    friend std::ostream& operator<<(std::ostream& out, gr&);
};


#endif /* yagbv_hpp */
