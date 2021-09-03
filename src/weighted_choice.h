#ifndef DENOVONEAR_WEIGHTED_CHOICE_H_
#define DENOVONEAR_WEIGHTED_CHOICE_H_

#include <random>
#include <vector>
#include <string>

struct AlleleChoice {
    int pos;
    char ref;
    char alt;
    double prob;
    
    // allow for splice altering sites outside the CDS
    int offset;
};

class Chooser {
    std::vector<AlleleChoice> sites;
    std::vector<double> cumulative;
    std::uniform_real_distribution<double> dist;
    std::mt19937_64 generator;
    void reset_sampler();

 public:
    Chooser();
    void add_choice(int site, double prob, char ref='N', char alt='N', int offset=0);
    int sampled_index();
    int choice_pos_only();
    AlleleChoice choice();
    double get_summed_rate();
    int len() { return sites.size() ;};
    AlleleChoice iter(int pos) { return sites[pos]; };
    void append(Chooser other);
};

#endif  // DENOVONEAR_WEIGHTED_CHOICE_H_
