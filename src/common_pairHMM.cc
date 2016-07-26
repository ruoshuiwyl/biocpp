//
// Created by ruoshui on 7/24/16.
//

#include "common_pairHMM.h"
#include "util.h"

double CommonPairHMM::ComputeReadLikelihood(std::vector<int8_t> &haplotype_bases, std::vector<int8_t> &read_bases,
                                            std::vector<int8_t> &read_qual, std::vector<int8_t> &read_insert_qual,
                                            std::vector<int8_t> &read_delete_qual, std::vector<int8_t> &overall_gcp,
                                            int32_t haplotype_index) {
    InitializeProbabilities(read_insert_qual, read_delete_qual, overall_gcp);
    InitializePriors(haplotype_bases,read_bases,read_qual,haplotype_index);
    double init_delete_value = 1.0 / haplotype_bases.size();
    for( int j = 0; j < pad_haplotype_length_; ++j){
        delete_matrix_[0][j] = init_delete_value;
    }
    for( int i = 1; i < pad_read_length_; ++i ) {
        for (int j = haplotype_index + 1; j < pad_haplotype_length_; ++j) {
            match_matrix_[i][j] = prior_[i][j] * ( match_matrix_[i-1][j-1] * transition_[i][MATCH2MATCH] +
                    insert_matrix_[i-1][j-1] * transition_[i][INSERT2MATCH] + delete_matrix_[i-1][j-1] * transition_[i][DELETE2MATCH] );
            insert_matrix_[i][j] = match_matrix_[i-1][j] * transition_[i][MATCH2INSERT] + insert_matrix_[i-1][j] * transition_[i][INSERT2INSERT];
            delete_matrix_[i][j] = match_matrix_[i][j-1] * transition_[i][MATCH2DELETE] + delete_matrix_[i][j-1] * transition_[i][DELETE2DELETE];
        }
    }
    double result = 0.0;
    for (int j = 1; j < pad_haplotype_length_; ++j) {
        result += match_matrix_[pad_read_length_-1][j] + insert_matrix_[pad_read_length_][j];
    }
    return result;
}

void CommonPairHMM::InitializeProbabilities(std::vector<int8_t> &read_insert_qual,
                                            std::vector<int8_t> &read_delete_qual,
                                            std::vector<int8_t> &overall_gcp) {
    int read_length = read_insert_qual.size();
    for( int i = 0; i != read_length; ++i ){
        transition_[i+1][MATCH2MATCH]  = 1.0 - (Util::ProbError(read_insert_qual[i]) - Util::ProbError(read_delete_qual[i]));
        transition_[i+1][MATCH2INSERT] = Util::ProbError(read_insert_qual[i]);
        transition_[i+1][MATCH2DELETE] = Util::ProbError(read_delete_qual[i]);
        transition_[i+1][INSERT2MATCH] = 1.0 - Util::ProbError(overall_gcp[i]);
        transition_[i+1][DELETE2MATCH] = 1.0 - Util::ProbError(overall_gcp[i]);
        transition_[i+1][INSERT2INSERT]= Util::ProbError(overall_gcp[i]);
        transition_[i+1][DELETE2DELETE]= Util::ProbError(overall_gcp[i]);
    }
}

void CommonPairHMM::InitializePriors(std::vector<int8_t> &haplotype_bases, std::vector<int8_t> &read_bases,
                                     std::vector<int8_t> &read_qual, int haplotype_index) {
    for( int i = 0; i < read_bases.size(); ++i ){
        const int8_t x = read_bases[i];
        const int8_t q = read_qual[i];
        for( int j = haplotype_index; j != haplotype_bases.size(); ++j) {
            const int8_t  y = haplotype_bases[j];
            prior_[i+1][j+1] = ( x == y || x == 'N' || y == 'N' ) ?
                               1.0 - Util::Prob(q) : Util::ProbError(q) / 3;
        }
    }
}