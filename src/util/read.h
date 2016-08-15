//
// Created by ruoshui on 7/24/16.
//

#ifndef BIOCPP_READ_H
#define BIOCPP_READ_H


#include <vector>
#include <cstdint>

class Read {
public:
    Read(){};
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
    Read(const std::string &bases, const std::string &read_quals, const std::string &read_ins_quals, const std::string &read_del_quals, const std::string read_gcp_quals) {
        for (int i = 0; i < bases.size(); ++i){
            bases_.push_back(bases[i]);
            base_qualities_.push_back(read_quals[i]-kQualityOffset);
            insert_qualities_.push_back(read_ins_quals[i]-kQualityOffset);
            delete_qualities_.push_back(read_del_quals[i]-kQualityOffset);
            gcp_qualities_.push_back(read_gcp_quals[i]-kQualityOffset);
        }
    }
    int32_t GetLength() {
        return bases_.size();
    }

    void set_bases(std::string &bases){
        for( int i = 0; i < bases.size(); ++i ) {
            bases_.push_back(bases[i]);
        }
    }
    void set_quals(std::string &quals) {
        for( int i = 0; i < quals.size(); ++i ) {
            base_qualities_.push_back(quals[i] - kQualityOffset);
        }
    }
    void set_ins_quals(std::vector<char> &ins_quals){
        insert_qualities_ = ins_quals;
    }
    void set_del_quals(std::vector<char> &del_quals){
        delete_qualities_ = del_quals;
    }
    void set_gcp_quals(std::vector<char> &gcp_qualitis) {
        gcp_qualities_ = gcp_qualitis;
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
    const int kQualityOffset = 33;
};


//class Haplotype{
//public:
//    Haplotype(const std::string &haplotype) {
//        haplotype_ = haplotype;
//    }
//    Haplotype(const std::vector<char> &haplotype){
//        haplotype_ = std::string(haplotype.begin(), haplotype.end());
//    }
//    const int size(){
//        haplotype_.size();
//    }
//
//private:
//    std::string haplotype_;
//
//};

#endif //BIOCPP_READ_H
