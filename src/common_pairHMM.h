//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_COMMON_PAIRHMM_H
#define BIOCPP_COMMON_PAIRHMM_H


#include "pairHMM.h"

class CommonPairHMM : public PairHMM {

private:
     double ComputeReadLikelihood(std::vector<int8_t> &haplotype_bases, std::vector<int8_t> &read_bases,
                                       std::vector<int8_t> &read_qual, std::vector<int8_t> &read_insert_qual,
                                       std::vector<int8_t> &read_delete_qual, std::vector<int8_t> &overall_gcp,
                                       int32_t haplotype_index);

    void InitializeProbabilities(std::vector<int8_t> &read_insert_qual,
                                 std::vector<int8_t> &read_delete_qual,
                                 std::vector<int8_t> &overall_gcp);

    virtual void InitializePriors(std::vector<int8_t> &haplotype_bases,
                                  std::vector<int8_t> &read_bases,
                                  std::vector<int8_t> &read_qual,
                                  int haplotype_index);

};

#endif //BIOCPP_COMMON_PAIRHMM_H
