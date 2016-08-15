//
// Created by ruoshui on 7/24/16.
//

#include "Log10_PairHMM.h"
#include "../util/util.h"

Log10_PairHMM::Log10_PairHMM(std::vector<std::vector<char>> &haplotypes, std::vector<Read> &reads): PairHMM(haplotypes, reads) {}


void Log10_PairHMM::Initialization() {
    PairHMM::Initialization();
    match_matrix_.resize(pad_max_read_length_ * pad_max_haplotype_length_, -DBL_MAX);
    insert_matrix_.resize(pad_max_read_length_ * pad_max_haplotype_length_, -DBL_MAX);
    delete_matrix_.resize(pad_max_read_length_ * pad_max_haplotype_length_, -DBL_MAX);
}
void Log10_PairHMM::QualToTransProb(const std::vector<char> &read_insert_quals,
                                    const std::vector<char> &read_delete_quals,
                                    const std::vector<char> &read_gcp_quals) {
    int read_length = read_insert_quals.size();
    for( int i = 0; i != read_length; ++i ){
//        transition_[kMatchToMatch][i+1]  = 1.0 - (Quality::ProbError(read_insert_quals[i]) + Quality::ProbError(read_delete_quals[i]));
        transition_[kMatchToMatch][i+1] = matchToMatchProbLog10(read_insert_quals[i], read_delete_quals[i]);
        transition_[kMatchToInsert][i+1] = Quality::ProbErrorLog10(read_insert_quals[i]);
        transition_[kMatchToDelete][i+1] = Quality::ProbErrorLog10(read_delete_quals[i]);
        transition_[kIndelToMatch][i+1] =  Quality::ProbLog10(read_gcp_quals[i]);
        transition_[kInsertToInsert][i+1]= Quality::ProbErrorLog10(read_gcp_quals[i]);
        transition_[kDeleteToDelete][i+1]= Quality::ProbErrorLog10(read_gcp_quals[i]);
    }

}

double Log10_PairHMM::matchToMatchProbLog10(const char ins_qual, const char dle_qual) {
     int minQual;
     int maxQual;
    if (ins_qual <= dle_qual) {
        minQual = ins_qual;
        maxQual = dle_qual;
    } else {
        minQual = dle_qual;
        maxQual = ins_qual;
    }
    return log1p (- std::min(1.0 ,pow(10, MathUtils::ApproximateLog10SumLog10(-0.1 * minQual, -0.1 * maxQual)))) * INV_LN10 ;
}



double Log10_PairHMM::subComputeReadLikelihoodGivenHaplotype(const std::vector<char> &haplotype_bases,
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
                                                   Quality::ProbLog10(qual) : Quality::ProbErrorLog10(qual);
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
            match_matrix_[GetMatrixIndex(i, j)] = prior_[GetMatrixIndex(i, j)]  + MathUtils::ApproximateLog10SumLog10 (
                    match_matrix_[GetMatrixIndex(i-1,j-1)]  + transition_[kMatchToMatch][i] ,
                    insert_matrix_[GetMatrixIndex(i-1, j-1)] +  transition_[kIndelToMatch][i] ,
                    delete_matrix_[GetMatrixIndex(i-1, j-1)] +  transition_[kIndelToMatch][i]
            );
            insert_matrix_[GetMatrixIndex(i, j)] = MathUtils::ApproximateLog10SumLog10(
                    match_matrix_[GetMatrixIndex(i-1, j)] + transition_[kMatchToInsert][i] ,
                    insert_matrix_[GetMatrixIndex(i-1,j)] + transition_[kInsertToInsert][i]
            );
            delete_matrix_[GetMatrixIndex(i, j)] = MathUtils::ApproximateLog10SumLog10(
                    match_matrix_[GetMatrixIndex(i, j-1)] + transition_[kMatchToDelete][i],
                    delete_matrix_[GetMatrixIndex(i, j-1)] + transition_[kDeleteToDelete][i]
            );
        }
    }

    const int32_t endI = pad_read_length_ - 1;
    double final_sum_prob = MathUtils::ApproximateLog10SumLog10(match_matrix_[GetMatrixIndex(endI, 0)], insert_matrix_[GetMatrixIndex(endI, 0)]);
    for( int i = 2; i < pad_haplotype_length_; ++i) {
        final_sum_prob += MathUtils::ApproximateLog10SumLog10( final_sum_prob, match_matrix_[GetMatrixIndex(endI, i)] ,  insert_matrix_[GetMatrixIndex(endI, i)]);
    }
    return final_sum_prob;
}


int Log10_PairHMM::ComputeLikeliHood(std::vector<double> &result) {
//    PairHMM::Initialization();
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
            std::vector<char> &haplotype_bases = haplotypes_[j];
            next_haplotype_index = prev_haplotype_bases_.empty() || prev_haplotype_bases_.size() !=haplotype_bases.size() ?
                                   0 : PairHMM::FoundNextHaplotypeIndex(prev_haplotype_bases_, haplotype_bases);
            double lk = subComputeReadLikelihoodGivenHaplotype(haplotype_bases, read_bases, read_quals, read_insert_quals,
                                                               read_delete_quals, read_gcp_quals, next_haplotype_index);
            result[GetMatrixIndex(i, j, haplotypes_.size())] = lk;
            prev_haplotype_bases_ = haplotypes_[j];
        }
    }
    return result.size();
}