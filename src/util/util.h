//
// Created by ruoshui on 7/25/16.
//

#ifndef BIOCPP_UTIL_H
#define BIOCPP_UTIL_H


#include <vector>
#include <cstdint>
#include <cmath>

 class Quality {
public:
    static void Initilization ();
    static double ProbError(int8_t qual);
    static double Prob(int8_t qual);
};

#endif //BIOCPP_UTIL_H
