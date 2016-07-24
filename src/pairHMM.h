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


protected:
    virtual double ComputeReadLikelihood(std::vector<int8_t> &haplotype_bases, std::vector<int8_t> &read_bases,
                                             std::vector<int8_t> &read_qual, std::vector<int8_t> &read_insert_qual,
                                             std::vector<int8_t> &read_delete_qual, std::vector<int8_t> &overall_gcp,
                                             int32_t haplotype_index);

private:
    std::vector<std::vector<int8_t>> haplotypes;
    std::vector<Read> reads_;
    std::vector<double> likelihood_;
    int32_t max_haplotype_length_;
    int32_t max_read_length_;
    int32_t pad_max_haplotype_length_;
    int32_t pad_max_read_length_;
    std::vector<std::vector<double>> match_matrix_;
    std::vector<std::vector<double>> insert_matrix_;
    std::vector<std::vector<double>> delete_matrix_;
    std::vector<std::vector<double>> transition_;
    std::vector<std::vector<double>> prior_;
};


#endif //BIOCPP_PAIRHMM_H
