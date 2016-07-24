//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_READ_H
#define BIOCPP_READ_H


#include <vector>
#include <cstdint>

class Read {
public:
    Read (std::vector<uint8_t> bases_,
          std::vector<uint8_t> base_qualities_,
          std::vector<uint8_t> insert_qualities_,
          std::vector<uint8_t> delete_qualities_,
          std::vector<uint8_t> gcp_qualities_) :
            bases_length_(bases_.size()),
            left_padding_(0),
            right_padding_(0),
            bases_(bases_),
            base_qualities_(base_qualities_),
            insert_qualities_(insert_qualities_),
            delete_qualities_(delete_qualities_),
            gcp_qualities_(gcp_qualities_) {

    };
private:
    int32_t bases_length_;
    int32_t left_padding_;
    int32_t right_padding_;
    std::vector<uint8_t> bases_;
    std::vector<uint8_t> base_qualities_;
    std::vector<uint8_t> insert_qualities_;
    std::vector<uint8_t> delete_qualities_;
    std::vector<uint8_t> gcp_qualities_;


};


#endif //BIOCPP_READ_H
