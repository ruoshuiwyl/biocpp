//
// Created by ruoshui on 7/25/16.
//

#ifndef BIOCPP_UTIL_H
#define BIOCPP_UTIL_H


#include <vector>

enum class CigarOperation {
    M = 0,  //Match or Mismatch
    I = 1,  //Insert
    D = 2,  //Delete
    N = 3,  //Skipped region from the reference.
    S = 4,  //Soft clip
    H = 5,  //Hard clip
    P = 6,  //padding
    EQ = 7, //Match
    X = 8,  //Mismatch
};

#include <vector>
#include <cstdint>
#include <cmath>
#include <cfloat>
#include <string>
#include <algorithm>
#include <bits/move.h>

class Quality {
public:
    static void Initilization ();
    static double ProbError(int8_t qual);
    static double Prob(int8_t qual);
     static double ProbErrorLog10(int8_t qual);
     static double ProbLog10(int8_t qual);
};

class MathUtils {
public:
    class JacobianLogTable {
    public:
        const static double MAX_TOLERANCE;
        static double get(const double difference) ;
        static void initialize() ;
        static int fastRound(double d) ;

    private:
        static double TABLE_STEP;
        static double INV_STEP;
        static std::vector<double> cache ;
    };
    static double ApproximateLog10SumLog10(std::vector<double> vals);
    static double ApproximateLog10SumLog10(double small, double big);
    static double ApproximateLog10SumLog10(double a1, double a2, double a3);

};
const std::string kCigarOpString = "MIDNSHP=X";

class CigarElem{
public:
    CigarElem(){}
    CigarElem(const int32_t len, const CigarOperation op): length_(len), operation_(op) {
    }
    CigarElem(const CigarElem &other): length_(other.length_), operation_(other.operation_){}
    CigarElem(const CigarElem &&other): length_(std::move(other.length_)), operation_(std::move(other.operation_)){}
    const CigarElem& operator=(const CigarElem &rhs) {
        length_ = rhs.length_;
        operation_ = rhs.operation_;
    }
    const int32_t length(){
        return length_;
    }
    const CigarOperation operation(){
        return operation_;
    }
    friend std::ostream& operator<<(std::ostream & os, const CigarElem &c);

private:
    int32_t length_;
    CigarOperation operation_;

};

typedef  std::vector<CigarElem> Cigar;

//const double kNeInf = -DBL_MAX;
#endif //BIOCPP_UTIL_H
