//
// Created by ruoshui on 7/24/16.
//

#include <cassert>
#include "pairHMM.h"
PairHMM::PairHMM(std::vector<std::vector<int8_t>> &haplotypes, std::vector<Read> &reads) :
        haplotypes_(haplotypes),
        reads_(reads) {
    assert( !haplotypes.empty() );
    assert( !reads_.empty());
    max_haplotype_length_ = haplotypes_[0].size();
    for (int i = 1; i !=  haplotypes_.size(); ++i) {
        if (haplotypes_[i].size() > max_haplotype_length_) {
            max_haplotype_length_ = haplotypes_[i].size();
        }
    }
    max_read_length_ = reads_[0].GetLength();
    for (int i = 1; i != reads_.size(); ++i) {
        if (reads[i].GetLength() > max_read_length_) {
            max_read_length_ = reads_[i].GetLength();
        }
    }
    pad_max_haplotype_length_ = max_haplotype_length_ + 1;
    pad_max_read_length_ = max_read_length_ + 1;
    Initialization();
}

void PairHMM::Initialization() {
    for ( int i = 0; i < pad_max_read_length_; ++i ) {
        std::vector<double> v(pad_max_haplotype_length_,0.0);
        match_matrix_.push_back(v);
        insert_matrix_.push_back(v);
        delete_matrix_.push_back(v);
        prior_.push_back(v);
    }
    for ( int i = 0; i < 7; ++i ) {
        transition_[i].reserve(pad_max_read_length_);
        for( int j = 0; j < pad_max_read_length_; ++j ) {
            transition_[i][j] = 0.0;
        }
    }
    for ( int i = 0; i < max_read_length_; ++i ){
        std::vector<double> v(max_haplotype_length_,0.0);
        likelihood_.push_back(v);
    }
}

std::vector<std::vector<double>> PairHMM::ComputeLikeliHood() {
    int32_t read_count = 0;
        for ( auto it = reads_.begin(); it != reads_.end(); ++it ){
            const std::vector<int8_t> &read_bases = it->GetReadBases();
            const std::vector<int8_t> &read_qual = it->GetReadBaseQualities();
            const std::vector<int8_t> &insert_qual = it->GetReadInsertQualities();
            const std::vector<int8_t> &delete_qual = it->GetReadDeleteQualities();
            const std::vector<int8_t> &gcp_qual = it->GetReadGCPQualities();
            bool is_first_haplotype = true;
            for( int j = 0; j != haplotypes_.size(); ++j ){
                std::vector<int8_t> *haplotype = &haplotypes_[j];
                std::vector<int8_t> *next_haplotype = j == haplotypes_.size() - 1 ? nullptr : &haplotypes_[j+1];
                double lk = ComputeReadLikelihood(haplotype, read_bases, read_qual, insert_qual, delete_qual, gcp_qual, is_first_haplotype, next_haplotype);
                likelihood_[read_count][j] = lk;
                is_first_haplotype = false;
            }
            read_count++;
        }
}
double PairHMM::ComputeReadLikelihood(const std::vector<int8_t> *haplotype_bases, const std::vector<int8_t> &read_bases,
                                      const std::vector<int8_t> &read_qual, const std::vector<int8_t> &read_insert_qual,
                                      const std::vector<int8_t> &read_delete_qual,
                                      const std::vector<int8_t> &overall_gcp, const bool recache_read_value,
                                      const std::vector<int8_t> *next_haplotype) {
    pad_haplotype_length_ = haplotype_bases->size() + 1;
    pad_read_length_ = read_bases.size();
    haplotype_index_ = recache_read_value ? 0 : haplotype_index_;
    int32_t next_haplotype_index = (next_haplotype == nullptr && haplotype_bases->size() != next_haplotype->size()) ? 0 : FindFristDiffHaplotype(haplotype_bases, next_haplotype);
    double result = SubComputeReadLikelihood(*haplotype_bases, read_bases, read_qual, read_insert_qual,
                                             read_delete_qual,
                                             overall_gcp, haplotype_index_, recache_read_value, next_haplotype_index);
    assert( result >= 0.0 );
    haplotype_index_ = haplotype_index_ < next_haplotype_index ? 0 : next_haplotype_index;
    return result;
}
int32_t PairHMM::FindFristDiffHaplotype(const std::vector<int8_t> *haploype,
                                        const std::vector<int8_t> *next_hapotype) {
    assert( haploype != nullptr );
    assert( next_hapotype != nullptr );
    assert( haploype->size() && next_hapotype->size());
    for( int i = 0; i < haploype->size() && i < next_hapotype->size(); ++i  ){
        if( haploype[i] != next_hapotype[i] ){
            return i;
        }
    }
    return haploype->size();
}