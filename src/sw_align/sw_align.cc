//
// Created by ruoshui on 8/17/16.
//

#include "sw_align.h"


Cigar SWAlignment::ComputeAlign(const std::string &reference, const std::string &alternate) {
    if( reference.size() + 1 > max_col_num_ || alternate.size() + 1> max_col_num_){
        max_col_num_ = ((std::max(reference.size(), alternate.size())) + 8 / 8) << 3;
        Initilization(max_col_num_);
    }
    int nrow = reference.size() + 1;
    int ncol = alternate.size() + 1;
    std::vector<int> best_gap_v;
    best_gap_v.resize(ncol+1, INT32_MIN / 2);
    std::vector<int> gap_size_v;
    gap_size_v.resize(ncol+1, 0);
    std::vector<int> best_gap_h;
    best_gap_h.resize(nrow+1, INT32_MIN / 2);
    std::vector<int> gap_size_h;
    gap_size_h.resize(nrow+1, 0);


    for (int i = 1; i < nrow; ++i ) {
        const char ref_base = reference[i-1];
        for (int j = 1; j < ncol; ++j) {
            const char alt_base = alternate[j-1];
            const int step_diag = WD(ref_base , alt_base) + sw_matrix_[GetMatrixIndex(i-1, j-1)];
            int prev_gap = sw_matrix_[GetMatrixIndex(i-1,j)] + kOpen;
            best_gap_v[j] += kExtend;
            if (prev_gap > best_gap_v[j]) {
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            } else {
                gap_size_v[j]++;
            }
            const int step_down = best_gap_v[j];
            const int kd = gap_size_v[j];
            prev_gap = sw_matrix_[GetMatrixIndex(i,j-1)] + kOpen;
            best_gap_h[i] += kExtend;
            if (prev_gap > best_gap_h[i]) {
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            } else {
                gap_size_h[i]++;
            }
            const int step_right = best_gap_h[i];
            const int ki = gap_size_h[i];
            if ( step_diag >= step_down && step_diag >= step_right) {
                sw_matrix_[GetMatrixIndex(i,j)] = step_diag;
                btrack_matrix_[GetMatrixIndex(i,j)] = 0;
            } else if ( step_right >=  step_down) {
                sw_matrix_[GetMatrixIndex(i,j)] = step_right;
                btrack_matrix_[GetMatrixIndex(i,j)] = -ki;

            } else  {
                sw_matrix_[GetMatrixIndex(i,j)] = step_down;
                btrack_matrix_[GetMatrixIndex(i,j)] = kd;
            }
//            std::cout << i << "\t" << j << "\t" <<  sw_matrix_[GetMatrixIndex(i, j)] << std::endl;
        }
    }
    return ComputeCigar( reference.size(), alternate.size());
}

Cigar SWAlignment::ComputeAlign(const std::vector<char> &reference, std::vector<char> &alternate) {
    std::string ref = reference.data();
    std::string alt = alternate.data();
    return  ComputeAlign(ref, alt);
}

void SWAlignment::Initilization(const int col_num) {
    sw_matrix_.resize(col_num * col_num, 0);
    btrack_matrix_.resize(col_num * col_num, 0);
}

const int SWAlignment::GetMatrixIndex(const int row, const int col) {
    return row * max_col_num_ + col;
}

Cigar SWAlignment::ComputeCigar(const int ref_len, const int alt_len) {
    CigarElem cigar_elem;
    int p1 = 0, p2 = alt_len, segment_length = 0;
    int max_score = INT32_MIN;

    for (int i = 1; i <= ref_len; ++i) {
        int cur_score = sw_matrix_[GetMatrixIndex(i, alt_len)];
        if (max_score < cur_score) {
            max_score = cur_score;
            p1 = i;
        }
    }

    for (int i = 1; i <= alt_len; ++i) {
        int cur_score = sw_matrix_[GetMatrixIndex(ref_len, i)];
        if (cur_score > max_score || (cur_score == max_score && abs(ref_len - i) < abs(p1 - p2))) {
            p1 = ref_len;
            p2 = i;
            max_score = cur_score;
            segment_length = alt_len - i;
        }
    }
    Cigar cigar;
    if (segment_length > 0) {
        cigar_elem = CigarElem(segment_length, CigarOperation::S);
        cigar.push_back(cigar_elem);
        segment_length  = 0;
    }
    CigarOperation state = CigarOperation::M;
    do {
        int btr = btrack_matrix_[GetMatrixIndex(p1,p2)];
        int step_length = 1;
        CigarOperation new_state;
        if ( btr > 0 ) {
            new_state = CigarOperation::D;
            step_length = btr;
        } else if ( btr < 0) {
            new_state = CigarOperation ::I;
            step_length = -btr;
        } else {
            new_state = CigarOperation::M;
        }
        switch ( new_state ) {
            case CigarOperation::D:{
                p1 -= step_length;
                break;
            }
            case CigarOperation::I: {
                p2 -= step_length;
                break;
            }
            case CigarOperation ::M: {
                p1--;
                p2--;
                break;
            }
        }
        if ( new_state == state ) {
            segment_length += step_length;
        } else {
            cigar_elem = CigarElem(segment_length,state);
            cigar.push_back(cigar_elem);
            segment_length = step_length;
            state = new_state;
        }

    } while (p1 > 0 && p2 > 0);
    cigar_elem = CigarElem(segment_length,state);
    cigar.push_back(cigar_elem);
    if ( p2 > 0 ) {
        cigar_elem = CigarElem(p2, CigarOperation::S);
        cigar.push_back(cigar_elem);
    }
    return cigar;
}