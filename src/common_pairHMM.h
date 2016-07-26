//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_COMMON_PAIRHMM_H
#define BIOCPP_COMMON_PAIRHMM_H


#include "pairHMM.h"

class CommonPairHMM : public PairHMM {

private:
    double SubComputeReadLikelihood(const std::vector<int8_t> &haplotype_bases, const std::vector<int8_t> &read_bases,
                                        const std::vector<int8_t> &read_qual, const std::vector<int8_t> &read_insert_qual,
                                        const std::vector<int8_t> &read_delete_qual, const std::vector<int8_t> &overall_gcp,
                                        const int32_t haplotype_index, const bool recacheReadValues,
                                        const int32_t next_haplotype_index);

    void InitializeProbabilities(const std::vector<int8_t> &read_insert_qual,
                                 const std::vector<int8_t> &read_delete_qual,
                                 const std::vector<int8_t> &overall_gcp);

    virtual void InitializePriors(const std::vector<int8_t> &haplotype_bases,
                                  const std::vector<int8_t> &read_bases,
                                  const std::vector<int8_t> &read_qual,
                                  const int32_t haplotype_index);

};

#endif //BIOCPP_COMMON_PAIRHMM_H
