//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_PAIRHMM_H
#define BIOCPP_PAIRHMM_H


#include <string>
#include <vector>
#include "read.h"



class PairHMM {
public:
    PairHMM( std::vector<std::vector<int8_t>> &haplotypes, std::vector<Read> &reads);
    enum Type {
        MATCH2MATCH = 0,
        MATCH2INSERT,
        MATCH2DELETE,
        INSERT2MATCH,
        DELETE2MATCH,
        INSERT2INSERT,
        DELETE2DELETE
    };


protected:
    virtual double ComputeReadLikelihood(std::vector<int8_t> &haplotype_bases, std::vector<int8_t> &read_bases,
                                             std::vector<int8_t> &read_qual, std::vector<int8_t> &read_insert_qual,
                                             std::vector<int8_t> &read_delete_qual, std::vector<int8_t> &overall_gcp,
                                             int32_t haplotype_index);
    virtual void InitializeProbabilities(std::vector<int8_t> &read_insert_qual, std::vector<int8_t> &read_delete_qual, std::vector<int8_t> &overall_gcp);
    virtual void InitializePriors(std::vector<int8_t> &haplotype_bases, std::vector<int8_t> &read_bases, std::vector<int8_t> &read_qual, int haplotype_index)

protected:
    std::vector<std::vector<int8_t>> &haplotypes_;
    std::vector<Read> reads_;
    std::vector<double> likelihood_;
    int32_t max_haplotype_length_;
    int32_t max_read_length_;
    int32_t pad_max_haplotype_length_;
    int32_t pad_max_read_length_;
    int32_t pad_read_length_;
    int32_t pad_haplotype_length_;
    std::vector<std::vector<double>> match_matrix_;
    std::vector<std::vector<double>> insert_matrix_;
    std::vector<std::vector<double>> delete_matrix_;
    std::vector<double> transition_[7];
    std::vector<std::vector<double>> prior_;
};


#endif //BIOCPP_PAIRHMM_H
