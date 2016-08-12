//
// Created by ruoshui on 7/24/16.
//

#include "common_pairHMM.h"
#include "../util/util.h"





int CommonPairHMM::ComputeLikeliHood(std::vector<double> &result) {
    Initialization();
    result.resize(reads_.size() * haplotypes_.size(), 0.0);

    for ( int i = 0; i < reads_.size(); ++i ) {
        std::vector<char> &read_bases = reads_[i].GetReadBases();
        std::vector<char> &read_quals = reads_[i].GetReadBaseQualities();
        std::vector<char> &read_insert_quals = reads_[i].GetReadInsertQualities();
        std::vector<char> &read_delete_quals = reads_[i].GetReadDeleteQualities();
        std::vector<char> &read_gcp_quals = reads_[i].GetReadGCPQualities();
        PairHMM::prev_haplotype_bases_.clear();
        QualToTransProb(read_insert_quals, read_delete_quals, read_gcp_quals);
        int next_haplotype_index = 0;
        for( int j = 0; j < haplotypes_.size(); ++j){
            std::vector<char> &haplotype_bases = haplotypes_[i];
            next_haplotype_index = PairHMM::FoundNextHaplotypeIndex(prev_haplotype_bases_,haplotype_bases, next_haplotype_index);
            double lk = subComputeReadLikelihoodGivenHaplotype(haplotype_bases, read_bases, read_quals, read_insert_quals,
                                                               read_delete_quals, read_gcp_quals, next_haplotype_index);
            result[GetMatrixIndex(i, j, reads_.size())] = lk;
            prev_haplotype_bases_ = haplotypes_[i];
        }
    }
    return result.size();
}

void CommonPairHMM::QualToTransProb(const std::vector<char> &read_insert_quals,
                                    const std::vector<char> &read_delete_quals,
                                    const std::vector<char> &read_gcp_quals) {
    int read_length = read_insert_quals.size();
    for( int i = 0; i != read_length; ++i ){
        transition_[i+1][kMatchToMatch]  = 1.0 - (Quality::ProbError(read_insert_quals[i]) + Quality::ProbError(read_delete_quals[i]));
        transition_[i+1][kMatchToInsert] = Quality::ProbError(read_insert_quals[i]);
        transition_[i+1][kMatchToDelete] = Quality::ProbError(read_delete_quals[i]);
        transition_[i+1][kIndelToMatch] = 1.0 - Quality::ProbError(read_gcp_quals[i]);
        transition_[i+1][kInsertToInsert]= Quality::ProbError(read_gcp_quals[i]);
        transition_[i+1][kDeleteToDelete]= Quality::ProbError(read_gcp_quals[i]);
    }
}

double CommonPairHMM::subComputeReadLikelihoodGivenHaplotype(const std::vector<char> &haplotype_bases,
                                                             const std::vector<char> &read_bases,
                                                             const std::vector<char> &read_quals,
                                                             const std::vector<char> &read_insert_quals,
                                                             const std::vector<char> &read_delete_quals,
                                                             const std::vector<char> &read_gcp_quals,
                                                             const int32_t hap_start_index) {
    for (int i = 0; i < read_bases.size(); ++i) {
        const char rbase = read_bases[i];
        const char qual = read_quals[i];
        for (int j = 0; j < haplotype_bases.size(); ++j) {
            const char hbase = haplotype_bases[j];
            prior_[GetMatrixIndex(i + 1, j + 1)] = ( rbase == hbase || rbase == 'N' || hbase == 'N' ) ?
                                                   Quality::Prob(qual) : Quality::ProbError(qual);
        }
    }
    pad_haplotype_length_ = haplotype_bases.size() + 1;
    if (prev_haplotype_bases_.empty() || prev_haplotype_bases_.size() != haplotype_bases.size()) {
        const double delete_initial_value = log10(1.0/ haplotype_bases.size());
        for (int i = 0; i < pad_haplotype_length_; ++i) {
            delete_matrix_[GetMatrixIndex(0, i)] = delete_initial_value;
        }
    }
    pad_read_length_ = read_bases.size() + 1;



    for (int i = 1; i < pad_read_length_; ++i) {
        for (int j = hap_start_index + 1; j < pad_haplotype_length_; ++j) {
            match_matrix_[GetMatrixIndex(i, j)] = prior_[GetMatrixIndex(i, j)] * (
                    match_matrix_[GetMatrixIndex(i-1,j-1)] * transition_[kMatchToMatch][i] +
                    insert_matrix_[GetMatrixIndex(i-1, j-1)] * transition_[kIndelToMatch][i] +
                    delete_matrix_[GetMatrixIndex(i-1, j-1)] * transition_[kIndelToMatch][i]
            );
            insert_matrix_[GetMatrixIndex(i, j)] = (
                    match_matrix_[GetMatrixIndex(i-1, j)] * transition_[kMatchToInsert][i] +
                    insert_matrix_[GetMatrixIndex(i-1,j)] * transition_[kInsertToInsert][i]
            );
            delete_matrix_[GetMatrixIndex(i, j)] = (
                    match_matrix_[GetMatrixIndex(i, j-1)] * transition_[kMatchToDelete][i],
                    delete_matrix_[GetMatrixIndex(i, j-1)] * transition_[kDeleteToDelete][i]
            );
        }
    }

    const int32_t endI = pad_read_length_ - 1;
    double final_sum_prob = 0.0 ;
    for( int i = 1; i < pad_haplotype_length_; ++i) {
        final_sum_prob += (match_matrix_[GetMatrixIndex(endI, i)] + insert_matrix_[GetMatrixIndex(endI, i)]);
    }
    return final_sum_prob;

}
