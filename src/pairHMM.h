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
    void Initialization();
    std::vector<std::vector<double>>  ComputeLikeliHood();

protected:
    double ComputeReadLikelihood(const std::vector<int8_t> *haplotype_bases, const std::vector<int8_t> &read_bases,
                                         const std::vector<int8_t> &read_qual, const std::vector<int8_t> &read_insert_qual,
                                         const std::vector<int8_t> &read_delete_qual, const std::vector<int8_t> &overall_gcp,
                                         const bool recache_read_value, const std::vector<int8_t> *next_haplotype);
    int32_t FindFristDiffHaplotype(const std::vector<int8_t> *haploype, const std::vector<int8_t> *next_hapotype);
    virtual double SubComputeReadLikelihood(const std::vector<int8_t> &haplotype_bases, const std::vector<int8_t> &read_bases,
                                                const std::vector<int8_t> &read_qual, const std::vector<int8_t> &read_insert_qual,
                                                const std::vector<int8_t> &read_delete_qual, const std::vector<int8_t> &overall_gcp,
                                                const int32_t haplotype_index, const bool recacheReadValues,
                                                const int32_t next_haplotype_index);

    virtual void InitializeProbabilities(const std::vector<int8_t> &read_insert_qual,
                                         const std::vector<int8_t> &read_delete_qual,
                                         const std::vector<int8_t> &overall_gcp);
    virtual void InitializePriors(const std::vector<int8_t> &haplotype_bases, const std::vector<int8_t> &read_bases,
                                  const std::vector<int8_t> &read_qual, const int32_t haplotype_index);

protected:
    std::vector<std::vector<int8_t>> &haplotypes_;
    std::vector<Read> reads_;
    std::vector<std::vector<double>> likelihood_;
    int32_t max_haplotype_length_;
    int32_t max_read_length_;
    int32_t pad_max_haplotype_length_;
    int32_t pad_max_read_length_;
    int32_t pad_read_length_;
    int32_t pad_haplotype_length_;
    int32_t haplotype_index_;
    std::vector<std::vector<double>> match_matrix_;
    std::vector<std::vector<double>> insert_matrix_;
    std::vector<std::vector<double>> delete_matrix_;
    std::vector<double> transition_[7];
    std::vector<std::vector<double>> prior_;
};


#endif //BIOCPP_PAIRHMM_H
