//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_PAIRHMM_H
#define BIOCPP_PAIRHMM_H


#include <string>
#include <vector>
#include "../util/read.h"


class PairHMM {
public:
    PairHMM( std::vector<std::vector<char>> &haplotypes, std::vector<Read> &reads);
    virtual int ComputeLikeliHood(std::vector<double> &result) = 0;


protected:
    virtual  void Initialization();
    virtual  void QualToTransProb(const std::vector<char> &read_insert_quals,
                                  const std::vector<char> &read_delete_quals,
                                  const std::vector<char> &read_gcp_quals) = 0;
    const int32_t FoundNextHaplotypeIndex(std::vector<char> &prev_hc, std::vector<char> &curr_hc, int curr_hc_index);
//    int32_t FindFristDiffHaplotype(const std::vector<char> *haploype, const std::vector<char> *next_hapotype);
    inline const int32_t GetMatrixIndex(int row , int col){
        return row * pad_max_haplotype_length_ + col;
    };
    inline const int32_t GetMatrixIndex(int row, int col, int col_num){
        return row * col_num + col;
    }
    std::vector<std::vector<char>> &haplotypes_;
    std::vector<Read> reads_;
    std::vector<std::vector<double>> likelihood_;
    int32_t max_haplotype_length_;
    int32_t max_read_length_;
    int32_t pad_max_haplotype_length_;
    int32_t pad_max_read_length_;
    int32_t pad_read_length_;
    int32_t pad_haplotype_length_;
//    int32_t haplotype_index_;
    const int32_t kPadSize = 1;
    const int32_t kMatchToMatch = 0;
    const int32_t kIndelToMatch = 1;
    const int32_t kMatchToInsert = 2;
    const int32_t kInsertToInsert = 3;
    const int32_t kMatchToDelete = 4;
    const int32_t kDeleteToDelete = 5;
    const static int32_t kTransitionMartixLength = 6;

    std::vector<double> match_matrix_;
    std::vector<double> insert_matrix_;
    std::vector<double> delete_matrix_;
    std::vector<double> transition_[kTransitionMartixLength];
    std::vector<double> prior_;
    
    std::vector<char> prev_haplotype_bases_;
};


#endif //BIOCPP_PAIRHMM_H
