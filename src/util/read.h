//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_READ_H
#define BIOCPP_READ_H


#include <vector>
#include <cstdint>

class Read {
public:
    Read (std::vector<char> &bases_,
          std::vector<char> &base_qualities_,
          std::vector<char> &insert_qualities_,
          std::vector<char> &delete_qualities_,
          std::vector<char> &gcp_qualities_) :
            bases_(bases_),
            base_qualities_(base_qualities_),
            insert_qualities_(insert_qualities_),
            delete_qualities_(delete_qualities_),
            gcp_qualities_(gcp_qualities_) {

    };
    int32_t GetLength() {
        return bases_.size();
    }
    std::vector<char> & GetReadBases() { return bases_;};
    std::vector<char> & GetReadBaseQualities() { return base_qualities_;};
    std::vector<char> & GetReadInsertQualities() { return insert_qualities_;};
    std::vector<char> & GetReadDeleteQualities() { return delete_qualities_;};
    std::vector<char> & GetReadGCPQualities() { return gcp_qualities_;};
private:
    std::vector<char> bases_;
    std::vector<char> base_qualities_;
    std::vector<char> insert_qualities_;
    std::vector<char> delete_qualities_;
    std::vector<char> gcp_qualities_;
};


#endif //BIOCPP_READ_H
