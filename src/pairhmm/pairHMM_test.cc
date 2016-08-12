//
// Created by ruoshui on 8/12/16.
//

#include <gtest/gtest.h>
#include <fstream>


static void TestCaseSetUp(const char *filename){
    std::ifstream in( filename);
    std::cout << in.is_open() << std::endl;
    int read_num, haplotype_num;
    in >> read_num >> haplotype_num;
    std::string read_bases, read_quals, read_ins_quals,  read_del_quals, read_gcp_quals;
    for (int i = 0; i < read_num; ++i ){
        in >> read_bases >> read_quals >> read_ins_quals >> read_del_quals >> read_gcp_quals;
        std::cout << read_bases << read_quals << read_ins_quals << read_del_quals << read_gcp_quals;
    }
    std::string halpotype;
    for (int i = 0; i < haplotype_num; ++i) {
        in >> halpotype;
    }
}

TEST(PairHMM, CommonPairHMM){
    TestCaseSetUp("data/tiny.in");




}