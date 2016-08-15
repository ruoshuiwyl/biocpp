//
// Created by ruoshui on 7/25/16.
//

#include "util.h"

static std::vector<double> prob_error_;
static std::vector<double> prob_;
static std::vector<double> prob_error_log10_;
static std::vector<double> prob_log10_;
static bool init = false;

void Quality::Initilization() {
    if (init == false) {
        prob_error_.reserve(128);
        prob_.reserve(128);
        prob_log10_.reserve(128);
        prob_error_log10_.reserve(128);
        for (int i = 0; i != 128; ++i) {
            prob_error_[i] = pow(10.0, -i / 10.0);
            prob_[i] = 1.0 - prob_error_[i];
            prob_error_log10_[i] = -0.1 * i;
            prob_log10_[i] = log10(1.0 - prob_error_[i]);
        }
        init = true;
    }
}

double Quality::Prob(int8_t qual) {
    Initilization();
    return prob_[qual];
}

double Quality::ProbError(int8_t qual) {
    Initilization();
    return prob_error_[qual];
}

double Quality::ProbErrorLog10(int8_t qual) {
    Initilization();
    return prob_error_log10_[qual];
}

double Quality::ProbLog10(int8_t qual) {
    Initilization();
    return prob_log10_[qual];
}

double MathUtils::JacobianLogTable::TABLE_STEP = 0.0001;
double MathUtils::JacobianLogTable::INV_STEP = 1.0 / TABLE_STEP;
const double MathUtils::JacobianLogTable::MAX_TOLERANCE = 8.0;

std::vector<double> MathUtils::JacobianLogTable::cache;

double MathUtils::JacobianLogTable::get(const double difference) {
    {
        if (cache.empty())
            initialize();
        int index = fastRound(difference * INV_STEP);
        return cache[index];
    }

}

void MathUtils::JacobianLogTable::initialize() {
        if (cache.empty()) {
            const int tableSize = (int) (MAX_TOLERANCE / TABLE_STEP) + 1;
            cache.resize(tableSize);
            for (int k = 0; k < tableSize; k++)
                cache[k] = log10(1.0 + pow(10.0, -((double) k) * TABLE_STEP));
        }
}

int MathUtils::JacobianLogTable::fastRound(double d) {
    return (d > 0.0)?  int(d + 0.5) : int( d - 0.5);
}


double MathUtils::ApproximateLog10SumLog10(double a1, double a2, double a3) {
    return ApproximateLog10SumLog10(a1, ApproximateLog10SumLog10(a2, a3));
}

double MathUtils::ApproximateLog10SumLog10(double small, double big) {
    if (small > big) {
        const double t = big;
        big = small;
        small = t;
    }

    if (small == -DBL_MAX || big == -DBL_MAX)
        return big;

    const double diff = fabs(big - small);
    if (diff >= JacobianLogTable::MAX_TOLERANCE)
        return big;
    return big + JacobianLogTable::get(diff);
}

double MathUtils::ApproximateLog10SumLog10(std::vector<double> vals) {
    double approx_sum = vals[0];
    int max_index = 0;
    for (int i = 1; i < vals.size(); ++i ){
        if (approx_sum < vals[i] ){
            approx_sum = vals[i];
            max_index = i;
        }
    }
    for (int i = 0; i < vals.size(); ++i) {
        if (i == max_index || vals[i] == DBL_MIN ) {
            continue;
        }
        const double diff = approx_sum - vals[i];
        if ( diff < JacobianLogTable::MAX_TOLERANCE ) {
            approx_sum += JacobianLogTable::get(diff);
        }

    }
    return approx_sum;
}