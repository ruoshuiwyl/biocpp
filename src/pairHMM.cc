//
// Created by ruoshui on 7/24/16.
//

#include "pairHMM.h"
PairHMM::PairHMM(std::vector<std::vector<int8_t>> &haplotypes, std::vector<Read> &reads) :
        haplotypes_(haplotypes),
        reads_(reads) {

}