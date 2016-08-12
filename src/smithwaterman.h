//
// Created by ruoshui on 7/2/16.
//

#ifndef BIOCPP_SMITHWATERMAN_H
#define BIOCPP_SMITHWATERMAN_H


#include <cstdint>
#include <string>
#include <vector>

struct ALignment{
    int32_t sw_score;
    int32_t sw_score_next_score;
    int32_t ref_begin;
    int32_t ref_end;
    int32_t query_begin;
    int32_t query_end;
    std::string cigar_str;
    std::vector<uint32_t> cigar;
};



class SWAlignment {
public:
    SWAlignment(){};
    ~SWAlignment(){};
    int32_t SWAlign(const std::string &ref, const std::string &query, ALignment &sw_align)  ;
    int32_t SWAlign(const char *ref, const int32_t rlen, const char *query, const int32_t qlen, ALignment &sw_align) ;

private:
    const int32_t kGapDeleteOpen = 0.0;
    const int32_t kGapInsertOpen = 0.0;
    const int32_t kGapExtension = 0.0;
    const int32_t kMatch = 0.0 ;
    const int32_t kMisMath = 0.0;
    int32_t  *score_matrix_;
    int32_t  score_matrix_size_;
    char  *ref_;
    char  *query_;
    int32_t ref_leng_;
    int32_t query_leng_;
    inline int32_t wd(const char x, const char y) const  {
        return x == y ? kMatch : kMisMath;
    }


};


#endif //BIOCPP_SMITHWATERMAN_H
