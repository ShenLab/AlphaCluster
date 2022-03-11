#ifndef DENOVONEAR_WEIGHTED_CHOICE_H_
#define DENOVONEAR_WEIGHTED_CHOICE_H_

#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <map>

struct AlleleChoice {
  int pos;
  char ref;
  char alt;
  double prob;
  // allow for splice altering sites outside the CDS
  int offset;
  double score;    
};

//class GeneralChooser{
//public:
//  virtual void generate_choices(int * positions_array, int positions_length) = 0;
//};

//class InheritedChooser : public GeneralChooser {
//std::uniform_int_distribution<int> dist;  
//  std::mt19937_64 generator;
//  int size;
//  std::vector<std::vector<int>> variants_by_individual;
//  void reset_sampler();
  
  
//public:
//  InheritedChooser();
//  void generate_choices(int * positions_array, int positions_length);  
//  void add_choice(std::vector<int> variants_for_individual);
//};

class Chooser {
    std::vector<AlleleChoice> sites;
    std::map<int,double> prob_per_site;
    std::vector<double> cumulative;
    std::vector<double> cumulative_score;  
    std::uniform_real_distribution<double> dist;
    std::uniform_int_distribution<int> dist_int;
    std::mt19937_64 generator;
    void reset_sampler();

    std::vector<std::vector<int>> variants_by_individual;
    std::vector<int> inherited_variants;
    int inherited_count;
    void reset_inherited_sampler();

public:
    Chooser();
    void generate_choices(int * positions_array, int positions_length);
    void generate_choices(int * positions_array, double * score_array, int positions_length);
    void generate_choices(std::vector<int> & positions_array,
			  std::vector<double> & score_array,
			  int positions_length);    
    void add_choice(int site, double prob, char ref='N', char alt='N', int offset=0);
    void add_choice(int site, double prob, char ref='N', char alt='N', int offset=0, double score = 1);  
    int sampled_index();
    int choice_pos_only();
    AlleleChoice choice();
    AlleleChoice choice(int chain);
    double get_summed_rate();
    double get_summed_score();
    double get_prob(int site){return prob_per_site[site];};
    std::map<int,double> get_prob_per_site(){return prob_per_site;};
    int len() { return sites.size() ;};
    std::vector<AlleleChoice> get_sites(){return sites;};
    AlleleChoice iter(int pos) { return sites[pos]; };
    void append(Chooser other);
    void generate_inherited_choices(int * positions_array, int positions_length);
    void add_inherited_choice(std::vector<int> variants_for_individual);
    void set_inherited_choices(std::vector<int> variants);
    std::vector<int> get_inherited_variants();
    int get_inherited_count() { return inherited_count;};
  //  void set_window(int start, int end);
};

#endif  // DENOVONEAR_WEIGHTED_CHOICE_H_
