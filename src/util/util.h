//
// Created by ruoshui on 7/25/16.
//

#ifndef BIOCPP_UTIL_H
#define BIOCPP_UTIL_H


#include <vector>
#include <cstdint>
#include <cmath>
#include <cfloat>

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

const double kNeInf = DBL_MIN;
#endif //BIOCPP_UTIL_H
