//
// Created by ruoshui on 8/12/16.
//

#include <gtest/gtest.h>
#include <fstream>
#include "../util/read.h"
#include "pairHMM.h"
#include "common_pairHMM.h"
#include "Log10_PairHMM.h"


std::vector<Read> reads;
std::vector<std::vector<char>> haplotypes;
static void TestCaseSetUp(const char *filename){
    std::ifstream in( filename);
    std::cout << in.is_open() << std::endl;
    int read_num, haplotype_num;
    in >> read_num >> haplotype_num;
    std::string read_bases, read_quals, read_ins_quals,  read_del_quals, read_gcp_quals;
    for (int i = 0; i < read_num; ++i ){
        in >> read_bases >> read_quals >> read_ins_quals >> read_del_quals >> read_gcp_quals;
//        std::cout << read_bases << read_quals << read_ins_quals << read_del_quals << read_gcp_quals;
        Read read(read_bases, read_quals, read_ins_quals, read_del_quals, read_gcp_quals);
        reads.push_back(read);
    }
    std::string halpotype;
    std::vector<char> haplotype_bases;
    for (int i = 0; i < haplotype_num; ++i) {
        in >> halpotype;
        haplotype_bases.assign(halpotype.begin(), halpotype.end());
        haplotypes.push_back(haplotype_bases);
    }
    in.close();
}

TEST(PairHMM, CommonPairHMM){
    TestCaseSetUp("/home/ruoshui/git/biocpp/src/data/tiny.in");
    std::vector<double> result;
//    PairHMM *pairhmm = new CommonPairHMM(haplotypes, reads);
    PairHMM *pairhmm = new Log10_PairHMM(haplotypes, reads);
    pairhmm->ComputeLikeliHood(result);
    for ( auto r : result) {
        std::cout << r << std::endl;
    }
    delete pairhmm;

}