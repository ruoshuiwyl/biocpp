//
// Created by ruoshui on 8/12/16.
//

#include <gtest/gtest.h>
#include <fstream>
#include "../util/read.h"
#include "pairHMM.h"
#include "common_pairHMM.h"
#include "Log10_PairHMM.h"


static std::vector<Read> reads;
static std::vector<std::vector<char>> haplotypes;
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
static void TestCaseSetup1(){
    const int read_num = 7;
    Read read[read_num];
    std::string read_bases[read_num] = {
            "GGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGG",
            "GGTACCCTCAGCCGGCCCGCTCGCCCGGGCCTGACCTGAGG",
            "CCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGG",
            "AGCCGGCCCGCTCGCCCGGGTCTGACCTGAGG",
            "AGCCGGCCCGCTCGCCCGGGTCTGACCTGAGG",
            "CCCGCCCGGGTCTGACCTGAGG",
            "CCGGGTCTGACCTGAGG"
    };
    std::string read_quals[read_num] = {
            ";@???555555555'5555555'555555'5555'555'55",
            ";A>>@@@B@@AA=9??A669A>648>'?7'<>A<5555555",
            "555555555355555555555555555555555555",
            "55555555555555555555555555555555",
            "55555555555555555555555555555555",
            ">>=8?<<8@@=>@A@@@AA@BA",
            ">?9>@>@AA@@@AB@C>"
    };
    std::vector<Read> reads;
    std::vector<char> ins_quals, del_quals, gcp_quals;
    for (int i = 0; i < read_num; ++i ){
        read[i].set_bases(read_bases[i]);
        read[i].set_quals(read_quals[i]);
        ins_quals.clear();
        del_quals.clear();
        gcp_quals.clear();
        for (int j = 0; j < read_bases[i].size() - 1; ++j) {
            ins_quals.push_back(40);
            del_quals.push_back(40);
            gcp_quals.push_back(10);
        }
        ins_quals.push_back(47);
        del_quals.push_back(45);
        gcp_quals.push_back(10);
        read[i].set_ins_quals(ins_quals);
        read[i].set_del_quals(del_quals);
        read[i].set_gcp_quals(gcp_quals);
        reads.push_back(read[i]);
    }
    std::vector<char> haplotype1, haplotype2;
    std::string h_bases[2] = {
            "GGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGG",
            "GGTACCCTCAGCCGGCCCGCTCGCCCGGGTCTGACCTGAGG"
    };
    haplotype1 = std::vector<char>(h_bases[0].begin(), h_bases[0].end());
    haplotype2 = std::vector<char>(h_bases[1].begin(), h_bases[1].end());

//    haplotype1.set_bases(h_bases[0]);
//    haplotype2.set_bases(h_bases[1]);
    std::vector<std::vector<char>> haplotypes;
//    reads.push_back(read);
    haplotypes.push_back(haplotype1);
    haplotypes.push_back(haplotype2);
}

TEST(PairHMM, CommonPairHMM){
//    TestCaseSetUp("/home/ruoshui/git/biocpp/src/data/tiny.in");
//    TestCaseSetup1();
    const int read_num = 7;
    Read read[read_num];
    std::string read_bases[read_num] = {
            "GGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGG",
            "GGTACCCTCAGCCGGCCCGCTCGCCCGGGCCTGACCTGAGG",
            "CCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGG",
            "AGCCGGCCCGCTCGCCCGGGTCTGACCTGAGG",
            "AGCCGGCCCGCTCGCCCGGGTCTGACCTGAGG",
            "CCCGCCCGGGTCTGACCTGAGG",
            "CCGGGTCTGACCTGAGG"
    };
    std::string read_quals[read_num] = {
            ";@???555555555'5555555'555555'5555'555'55",
            ";A>>@@@B@@AA=9??A669A>648>'?7'<>A<5555555",
            "555555555355555555555555555555555555",
            "55555555555555555555555555555555",
            "55555555555555555555555555555555",
            ">>=8?<<8@@=>@A@@@AA@BA",
            ">?9>@>@AA@@@AB@C>"
    };
    std::vector<Read> reads;
    std::vector<char> ins_quals, del_quals, gcp_quals;
    for (int i = 0; i < read_num; ++i ){
        read[i].set_bases(read_bases[i]);
        read[i].set_quals(read_quals[i]);
        ins_quals.clear();
        del_quals.clear();
        gcp_quals.clear();
        for (int j = 0; j < read_bases[i].size() - 1; ++j) {
            ins_quals.push_back(40);
            del_quals.push_back(40);
            gcp_quals.push_back(10);
        }
        ins_quals.push_back(47);
        del_quals.push_back(45);
        gcp_quals.push_back(10);
        read[i].set_ins_quals(ins_quals);
        read[i].set_del_quals(del_quals);
        read[i].set_gcp_quals(gcp_quals);
        reads.push_back(read[i]);
    }
    std::vector<char> haplotype1, haplotype2;
    std::string h_bases[2] = {
            "GGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGG",
            "GGTACCCTCAGCCGGCCCGCTCGCCCGGGTCTGACCTGAGG"
    };
    haplotype1 = std::vector<char>(h_bases[0].begin(), h_bases[0].end());
    haplotype2 = std::vector<char>(h_bases[1].begin(), h_bases[1].end());

//    haplotype1.set_bases(h_bases[0]);
//    haplotype2.set_bases(h_bases[1]);
    std::vector<std::vector<char>> haplotypes;
//    reads.push_back(read);
    haplotypes.push_back(haplotype1);
    haplotypes.push_back(haplotype2);
    std::vector<double> result;
    PairHMM *pairhmm = new CommonPairHMM(haplotypes, reads);
//    PairHMM *pairhmm = new Log10_PairHMM(haplotypes, reads);
    pairhmm->ComputeLikeliHood(result);
    for ( auto r : result) {
        std::cout << r << std::endl;
    }
    delete pairhmm;
}