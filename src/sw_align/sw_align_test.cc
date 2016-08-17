//
// Created by ruoshui on 8/17/16.
//

#include <gtest/gtest.h>
#include "sw_align.h"

TEST( GSWAlignment, sw){
    std::string ref = "ACACACTA";
    std::string alt = "AGCACACA";
    std::string ex_result = "1M1D5M1I1M";
    std::string ac_result;
    SWAlignment sw(2, -1, -1, -1);
    Cigar cigar =  sw.ComputeAlign(ref, alt);
    for( auto c : cigar) {

       ac_result.append(std::to_string(c.length()));
        ac_result.push_back(kCigarOpString[int( c.operation())]);
    }
    EXPECT_EQ(ex_result, ac_result);
    std::cout << std::endl;

}