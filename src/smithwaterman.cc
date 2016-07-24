//
// Created by ruoshui on 7/2/16.
//

#include "smithwaterman.h"
#include <algorithm>
#include <cstring>

int SWAlignment::SWAlign(const std::string &ref, const std::string &query, ALignment &sw_align)  {
    SWAlign(ref.c_str(), ref.size(), query.c_str(), query.size(), sw_align);
}
int SWAlignment::SWAlign(const char *ref,
                         const int32_t rlen,
                         const char *query,
                         const int32_t qlen,
                         ALignment &sw_align)  {

    if (rlen + 1 > score_matrix_size_ || qlen + 1 > score_matrix_size_) {
        score_matrix_size_ = rlen > qlen ? rlen + 1 : qlen + 1;
        delete[] score_matrix_;
        score_matrix_ = new int32_t[score_matrix_size_ * score_matrix_size_]{0};
    }

    int32_t *prev_row, *curr_row;
    int32_t best_gap_v[qlen+1]{INT32_MIN};
    int32_t best_size_v[qlen+1]{0};
    int32_t best_gap_h[rlen+1]{INT32_MIN};
    int32_t best_size_h[rlen+1]{0};

    //init

    for( int i = 1; i <= rlen; ++i ){
        prev_row = curr_row;
        curr_row = &score_matrix_[i * score_matrix_size_];
        const char rbase = ref[i-1];
        for (int j = 1; j <= qlen; ++j) {
            const char qbase = query[j-1];
            int32_t step_diag = prev_row[j-1] + wd(rbase, qbase);
            int32_t prev_gap = prev_row[j] + kGapDeleteOpen;
            best_gap_v[j] += kGapExtension;
            if (prev_gap > best_gap_v[j]) {
                best_gap_v[j] = prev_gap;
                best_size_v[j] = 1;
            } else {
                best_size_v[j]++;

            }
            int32_t step_down = best_gap_v[j];

            prev_gap = curr_row[j-1] + kGapInsertOpen;
            best_gap_h[i] += kGapExtension;
            if (prev_gap > best_gap_h[i]) {
                best_gap_h[i] = prev_gap;
                best_size_h[i] = 1;
            } else {
                best_size_h[i]++;
            }
            int32_t step_right = best_gap_h[i];
            bool diaghighest = (step_diag >= step_down) && step_diag >= step_right);
            if (diaghighest){
                curr_row[j] = step_diag > 0 ? step_diag : 0;
            } else if (step_right >= step_down) {
                curr_row[j] = step_right > 0 ? step_right : 0;
            } else {
                curr_row[j] = step_down > 0 ? step_down : 0;
            }
        }
    }

    //back trace
    int32_t max_score_index = rlen*score_matrix_size_+qlen;
    int32_t max_score = score_matrix_[rlen*score_matrix_size_+qlen];
    for (int i = 0; i < rlen; ++i) {
            if (score_matrix_[i*score_matrix_size_+qlen] > max_score) {
                max_score = score_matrix_[i*score_matrix_size_+qlen];
                max_score_index = i*score_matrix_size_;
            }
    }

    for( int begin = max_score_index, end = qlen; begin > 0; --begin)



}