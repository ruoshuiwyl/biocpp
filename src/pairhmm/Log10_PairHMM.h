//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_LOG10_PAIRHMM_H
#define BIOCPP_LOG10_PAIRHMM_H


#include "pairHMM.h"

class Log10_PairHMM : public PairHMM {
public:
    Log10_PairHMM(std::vector<std::vector<char>> &haplotypes, std::vector<Read> &reads);
    int ComputeLikeliHood(std::vector<double> &result);


private:
    void Initialization();
    void QualToTransProb(const std::vector<char> &read_insert_quals,
                         const std::vector<char> &read_delete_quals,
                         const std::vector<char> &read_gcp_quals);

    double matchToMatchProbLog10(const char ins_qual, const char dle_qual);

    double subComputeReadLikelihoodGivenHaplotype(const std::vector<char> &haplotype_bases,
                                                  const std::vector<char> &read_bases,
                                                  const std::vector<char> &read_quals,
                                                  const std::vector<char> &read_insert_quals,
                                                  const std::vector<char> &read_delete_quals,
                                                  const std::vector<char> &read_gcp_quals,
                                                  const int32_t hap_start_index);
    const double INV_LN10 = 0.43429448190325176;
};


#endif //BIOCPP_LOG10_PAIRHMM_H
