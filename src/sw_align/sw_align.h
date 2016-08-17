//
// Created by ruoshui on 8/17/16.
//

#ifndef  SRC_SW_ALIGN_H
#define  SRC_SW_ALIGN_H


#include <cstdint>
#include <vector>
#include <string>
#include "../util/util.h"

class SWAlignment {
public:
    SWAlignment(): kMatch(200), kMismatch(-150), kOpen(-260), kExtend(-11){
        Initilization(max_col_num_);

    };
    SWAlignment(const int match, const int mismatch, const int open, const int extend):
            kMatch(match), kMismatch(mismatch), kOpen(open), kExtend(extend){
        Initilization(max_col_num_);
    }
    Cigar ComputeAlign(const std::vector<char> &reference, std::vector<char> &alternate);
    Cigar ComputeAlign(const std::string &reference, const std::string &alternate);
private:
    void Initilization(const int col_num);
    const int GetMatrixIndex(const int row, const int col);
    inline int WD(const char ref_base, const char alt_base) {
        return ref_base == alt_base? kMatch : kMismatch;
    }
    Cigar ComputeCigar(const int ref_len, const int alt_len);
    const int kMatch;
    const int kMismatch;
    const int kOpen;
    const int kExtend;
    std::vector<int32_t> sw_matrix_;
    std::vector<int32_t> btrack_matrix_;
    int max_col_num_ = 500;

};



#endif
 //SRC_SW_ALIGN_H
