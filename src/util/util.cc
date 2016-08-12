//
// Created by ruoshui on 7/25/16.
//

#include "util.h"
static std::vector<double> prob_error_;
static std::vector<double> prob_;
static bool init = false;

void Quality::Initilization() {
    if ( init == false ) {
        prob_error_.reserve(128);
        prob_.reserve(128);
        for (int i = 0; i != 128; ++i) {
            prob_error_[i] = pow(10.0, -i / 10.0);
            prob_[i] = 1.0 - prob_error_[i];
        }
        init = true;
    }
}
double Quality::Prob(int8_t qual) {
    Initilization();
    return prob_[qual];
}

double Quality::ProbError(int8_t qual) {
    Initilization();
    return prob_error_[qual];
}