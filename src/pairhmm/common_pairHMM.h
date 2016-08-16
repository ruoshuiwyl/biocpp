//
// Created by ruoshui on 7/24/16.
//



#ifndef BIOCPP_COMMON_PAIRHMM_H
#define BIOCPP_COMMON_PAIRHMM_H


#include <cmath>
#include "pairHMM.h"

class CommonPairHMM : public PairHMM {
public:
    CommonPairHMM(std::vector<std::vector<char>> &haplotypes, std::vector<Read> &reads);

    int ComputeLikeliHood(std::vector<double> &result);

private:
    void Initialization();

    void QualToTransProb(const std::vector<char> &read_insert_quals,
                         const std::vector<char> &read_delete_quals,
                         const std::vector<char> &read_gcp_quals);

    double subComputeReadLikelihoodGivenHaplotype(const std::vector<char> &haplotype_bases,
                                                  const std::vector<char> &read_bases,
                                                  const std::vector<char> &read_quals,
                                                  const std::vector<char> &read_insert_quals,
                                                  const std::vector<char> &read_delete_quals,
                                                  const std::vector<char> &read_gcp_quals,
                                                  const int32_t hap_start_index);


    double kInitalConstant;
    double kLog10InitalConstant;
};

#endif //BIOCPP_COMMON_PAIRHMM_H
