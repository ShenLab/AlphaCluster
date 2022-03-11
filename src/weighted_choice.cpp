#include <random>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iterator>

#include "weighted_choice.h"
// Constructor for Chooser class - starts the random sampler
Chooser::Chooser() {
    /**
        Constructor for Chooser class
    */
    
    // start the random sampler
    std::random_device rd;
    generator.seed(rd());
    inherited_count = 0;
    inherited_variants = {};
}
 
// any time we add another choice, reset the sampler, so we can sample
// from all the possible entries.
void Chooser::reset_sampler() {
    std::uniform_real_distribution<double> temp(0.0, get_summed_rate());
    dist = temp;
}

// get index to randomly sampled element, using weighted probabilities
int Chooser::sampled_index() {
    if (cumulative.empty()) {
        return -1;
    }
    
    // get a random float between 0 and the cumulative sum
    double number = dist(generator);
    
    // figure out where in the list a random probability would fall
    auto pos = std::lower_bound(cumulative.begin(), cumulative.end(), number);
    return pos - cumulative.begin();
}

// adds another choice to the class object
//
// @param site site position (e.g. 100001)
// @param prob site mutation rate (e.g. 0.000000005)
// @param ref reference allele for site e.g. 'A'
// @param alt alternate allele for site e.g. 'T'
// @param offset number of bases the position is offset from the true chromosomal
//     site. This is relevant for non-CDS sites, such as within splice
//     regions.

void Chooser::add_choice(int site, double prob, char ref, char alt, int offset, double score) {
    // keep track of the cumulative sum for each added site
    double cumulative_sum = get_summed_rate() + prob;
    cumulative.push_back(cumulative_sum);

    double cumul_score = get_summed_score() + score;
    cumulative_score.push_back(cumul_score);
    
    sites.push_back(AlleleChoice {site, ref, alt, prob, offset, score});
    prob_per_site[site] += prob;
    reset_sampler();
}

void Chooser::generate_choices(int * positions_array, double * score_array, int positions_length){
  // int size = 0;
  AlleleChoice c;
  for(int i = 0; i < positions_length; i++){
    c = choice();
    positions_array[i] = c.pos;
    score_array[i] = c.score;
  }

}


void Chooser::generate_choices(std::vector<int>& positions,
			       std::vector<double>& scores,
			       int positions_length)
{
  AlleleChoice c;
  for(int i = 0; i < positions_length; i++){
    //select choice from correct chain
    c = choice();
    positions.push_back(c.pos);
    scores.push_back(c.score);
  }

}


void Chooser::generate_choices(int * positions_array, int positions_length){
  // int size = 0;
  for(int i = 0; i < positions_length; i++){
    positions_array[i] = choice().pos;
  }
}

// AlleleChoice Chooser::choice() {
//     /**
//         chooses a random element using a set of probability weights
        
//         @returns AlleleChoice struct containing the pos, ref and alt
//     */
    
//     if (cumulative.empty()) {
//       return AlleleChoice {-1, 'N', 'N', 0.0, 0,1};
//     }
    
//     // get a random float between 0 and the cumulative sum
//     double number = dist(generator);
    
//     // figure out where in the list a random probability would fall
//     auto pos = std::lower_bound(cumulative.begin(), cumulative.end(), number);
//     int offset = pos - cumulative.begin();
    
//     return sites[offset];
// }

// randomly sample a site, but only return the site position
// 
// This is a fast version for when we only want the site position, and want to
// avoid constructing a full AlleleChoice while returning.
int Chooser::choice_pos_only() {
    int idx = sampled_index();
    if (idx >= 0) {
        return sites[idx].pos;
    }
    return -1;
}

// randomly sample a site
AlleleChoice Chooser::choice() {
    int idx = sampled_index();
    if (idx >= 0) {
        return sites[idx];
    }
    return AlleleChoice {-1, 'N', 'N', 0.0, 0};
}

//  gets the cumulative sum for all the current choices.
double Chooser::get_summed_rate() {
    return (sites.empty()) ? 0.0 : cumulative.back() ;
}

//  gets the cumulative score for all the current choices.
double Chooser::get_summed_score() {
    return (sites.empty()) ? 0.0 : cumulative_score.back() ;
}


void Chooser::append(Chooser other) {
    
    double current = get_summed_rate();
    int len = other.sites.size();
    for (int i=0; i < len; i++) {
        cumulative.push_back(other.cumulative[i] + current);
        sites.push_back(other.sites[i]);
    }

    prob_per_site = other.get_prob_per_site();
    inherited_variants = other.get_inherited_variants();
    inherited_count = other.get_inherited_count();
    
    reset_sampler();
}

void Chooser::add_inherited_choice(std::vector<int> variants_for_individual){
  variants_by_individual.push_back(variants_for_individual);
  inherited_count += variants_for_individual.size();
  reset_sampler();
}

//void Chooser::set_window(int start, int end){}

void Chooser::reset_inherited_sampler() {
    // any time we add another choice, reset the sampler, so we can sample
    // from all the possible entries.
    std::uniform_int_distribution<int> temp(0, inherited_count);
    dist_int = temp;
}

void Chooser::set_inherited_choices(std::vector<int> variants) {
      inherited_variants = variants;
      inherited_count = inherited_variants.size();
}

std::vector<int> Chooser::get_inherited_variants() {
      return inherited_variants;
}

void Chooser::generate_inherited_choices(int * positions_array, int positions_length){
  std::vector<int> out; 
  //ForwardIterator out;
  std::sample(inherited_variants.begin(),
	      inherited_variants.end(),
	      std::back_inserter(out),
	      positions_length,
	      std::mt19937{std::random_device{}()});
  for(int i = 0; i < positions_length; i++){
    //    std::cout << i << std::endl;
    //std::cout<< i << " out = "<<  out[i] << std::endl;
    //std::cout<< i << " position = "<<  positions_array[i] << std::endl;
    positions_array[i] = out[i];
    //std::cout<< i << " position = "<<  positions_array[i] << std::endl;    
    
  }
  //std::copy(out.begin(), out.end(), positions_array);

  //for(int i = 0; i < positions_length; i++, out++){
  //  positions_array[i] = *out;
  //}
  
  /* int positions_idx = 0;
  int choice_idx = 0;
  
  int number;
  int * currentChoice;
  int currentChoiceSize;
  
  while (positions_idx < positions_length) {
    number = dist_int(generator);
    currentChoice = variants_by_individual[number].data();
    currentChoiceSize = variants_by_individual[number].size();
    while (positions_idx < positions_length && choice_idx < currentChoiceSize){
      positions_array[positions_idx] = currentChoice[choice_idx];
      positions_idx++;
      choice_idx++;
    }
  }
  */
}
