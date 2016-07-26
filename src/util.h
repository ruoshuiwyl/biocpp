//
// Created by ruoshui on 7/25/16.
//

#ifndef BIOCPP_UTIL_H
#define BIOCPP_UTIL_H


#include <vector>
#include <cstdint>
#include <cmath>

class Util {
public:
    static double ProbError(int8_t qual){
        if( init_ = false){
           prob_error_.reserve(128);
            prob_.reserve(128);
            for( int i = 0; i != 128; ++i ){
                prob_error_[i] = pow(10.0, - i / 10.0 );
                prob_[i] = 1 - prob_error_[i];
            }
           init_ = true;
        }
        double prob_error_[qual];
    }
    static double Prob(int8_t qual) {
        if( init_ = false){
            prob_error_.reserve(128);
            prob_.reserve(128);
            for( int i = 0; i != 128; ++i ){
                prob_error_[i] = pow(10.0, - i / 10.0 );
                prob_[i] = 1 - prob_error_[i];
            }
            init_ = true;
        }
        double prob_[qual];
    }
private:
    static std::vector<double> prob_error_;
    static std::vector<double> prob_;
    static bool init_ = false;
};


#endif //BIOCPP_UTIL_H
