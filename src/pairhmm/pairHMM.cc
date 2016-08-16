//
// Created by ruoshui on 7/24/16.
//

#include <cassert>
#include "pairHMM.h"
PairHMM::PairHMM(std::vector<std::vector<char>> &haplotypes, std::vector<Read> &reads) :
        haplotypes_(haplotypes),
        reads_(reads) {
    assert( !haplotypes.empty() );
    assert( !reads_.empty());
}

void PairHMM::Initialization() {
    max_haplotype_length_ = 0;
    for (int i = 0; i !=  haplotypes_.size(); ++i) {
        if (haplotypes_[i].size() > max_haplotype_length_) {
            max_haplotype_length_ = haplotypes_[i].size();
        }
    }
    max_read_length_ = 0;
    for (int i = 0; i != reads_.size(); ++i) {
        if (reads_[i].GetLength() > max_read_length_) {
            max_read_length_ = reads_[i].GetLength();
        }
    }
    pad_max_haplotype_length_ = max_haplotype_length_ + kPadSize;
    pad_max_read_length_ = max_read_length_ + kPadSize;

    prior_.resize(pad_max_haplotype_length_ * pad_max_read_length_, 0.0);
    for ( int i = 0; i < kTransitionMartixLength; ++i ) {
        transition_[i].resize(pad_max_read_length_, 0.0);
    }
}



const int32_t PairHMM::FoundNextHaplotypeIndex(std::vector<char> &prev_hc, std::vector<char> &curr_hc) {

    int index = 0;
    for ( int i = 0; i < curr_hc.size(); ++i ) {
        if (prev_hc[i] != curr_hc[i]) {
            index = i;
            break;
        }
    }

    return index;
}

