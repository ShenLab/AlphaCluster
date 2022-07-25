#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <assert.h>
#include <fstream>
#include <iterator>
#include "weighted_choice.h"

// gets the distances between all the pairs of elements from a list
//
// @param sites array of positions
// @return vector of paired positions
void _get_distances(int sites[],
		    double xyz[][3],
		    //const std::vector<std::vector<double>> & three_d_locations,
		    int & len,
		    double distances[]) {
    /**
        gets the distances between all the pairs of elements from a list
        
        @sites array of positions
        @return vector of paired positions
    */
    // get all non-repeating combinations of the sites
    int idx = 0;
    int codon_i, codon_j;
    
    for (int i=0; i < len; i++) {
      codon_i = floor((float)sites[i] / 3) + 1;
      for (int j=i+1; j < len; j++) {
	codon_j = floor((float)sites[j] / 3) + 1;
	distances[idx] = std::pow((xyz[codon_i][0] - xyz[codon_j][0])*(xyz[codon_i][0] - xyz[codon_j][0]) +
				  (xyz[codon_i][1] - xyz[codon_j][1])*(xyz[codon_i][1] - xyz[codon_j][1]) +
				  (xyz[codon_i][2] - xyz[codon_j][2])*(xyz[codon_i][2] - xyz[codon_j][2]), .5);
	
	idx += 1;
        }
    }
    if(distances[0] < 0){
      std::cout << distances[0] << std::endl;
      std::cout << sites[0] << std::endl;
      std::cout << sites[1] << std::endl;      
    }
}

void _get_distances(const std::vector<int> & sites,
		    const std::vector<std::vector<double>> & three_d_locations,
		    std::vector<double> & distances) {
    int size = sites.size();
    int len = ((size - 1) * size) / 2;
    int sites_array[size];
    for (int i=0; i<size; i++) {
        sites_array[i] = sites[i];
    }

    int aa_count = three_d_locations.size();
    double three_d_locations_array[aa_count][3];
    for (int i=0; i<aa_count; i++) {
        three_d_locations_array[i][0] = three_d_locations[i][0];
        three_d_locations_array[i][1] = three_d_locations[i][1];
        three_d_locations_array[i][2] = three_d_locations[i][2];	
    }

    double dist_array[len];
    _get_distances(sites_array, three_d_locations_array, size, dist_array);
    distances.resize(len);
    
    for (int i=0; i<len; i++) {
        distances[i] = dist_array[i];
    }
}
void _get_distances_coevol(const int sites[],
			   const double coevol[],
			   const int aa_len,
			   int & len,
			   double distances[])
{
  int idx = 0;
  int codon_i, codon_j;
  for (int i=0; i < len; i++) {
    codon_i = floor((float)sites[i] / 3) + 1;
    for (int j=i+1; j < len; j++) {
      codon_j = floor((float)sites[j] / 3);
      distances[idx] = coevol[codon_i*(aa_len-1) + codon_j];
      idx += 1;
    }
  }
}


void _get_distances_coevol(const std::vector<int> & sites,
			   const std::vector<double> & coevol,
			   std::vector<double> & distances)
{
    /**
       gets the distances between all the pairs of elements from a list
        
        @sites array of positions
        @return vector of paired positions
    */
    // get all non-repeating combinations of the sites
    int size = sites.size();
    int len = ((size - 1) * size) / 2;
    int sites_array[size];
    for (int i=0; i<size; i++) {
        sites_array[i] = sites[i];
    }
    double dist_array[len];
    
    _get_distances_coevol(sites_array,
			  coevol.data(),
			  (int) sqrt(coevol.size()),
			  size,
			  dist_array);    
 
    distances.resize(len);
    for (int i=0; i<len; i++) {
        distances[i] = dist_array[i];
    }
}


void _get_distances(int sites[], int & len, double distances[]) {
    /**
        gets the distances between all the pairs of elements from a list
        
        @sites array of positions
        @return vector of paired positions
    */
    // get all non-repeating combinations of the sites
    int idx = 0;
    for (int i=0; i < len; i++) {
        for (int j=i+1; j < len; j++) {
            // only include if the array positions differ, so we avoid finding
            // the distance to itself
            distances[idx] = abs(sites[i] - sites[j]);
            idx += 1;
        }
    }
}

void _get_distances(const std::vector<int> & sites, std::vector<double> & distances) {
    int size = sites.size();
    int len = ((size - 1) * size) / 2;
    int sites_array[size];
    for (int i=0; i<size; i++) {
        sites_array[i] = sites[i];
    }
    double dist_array[len];
    _get_distances(sites_array, size, dist_array);
    distances.resize(len);
    
    for (int i=0; i<len; i++) {
        distances[i] = dist_array[i];
    }
}

// check if any value in an vector is zero
//
// @param distances vector of values
// @return true/false for containing zero
bool _has_zero(double distances[], int & len) {
    /**
        @check if any value in an vector is zero
        
        @distances vector of values
        
        @return true/false for containing zero
    */
    return std::find(distances, distances + len, 0) != distances + len;
}

bool _has_zero(std::vector<double> distances) {
    double * dist = distances.data();
    int len = distances.size();
    return _has_zero(dist, len);
}

void _scale_distances(double distances[],
		      const double scores[],
		      int len)
{
  int idx = 0;
  for (int i=0; i < len; i++) {
    for (int j=i+1; j < len; j++) {
      distances[idx] = (1+distances[idx])/((scores[i]+scores[j])/2);
      idx += 1;
    }
  }
}

std::vector<double> _scale_distances(std::vector<double> distances,
		      const std::vector<double> scores)
{

  int size = distances.size();
  double dist_array[size];
  for (int i=0; i<size; i++) {
    dist_array[i] = distances[i];
  }
  
  _scale_distances(dist_array,
		   scores.data(),
		   scores.size());

  for (int i=0; i<size; i++) {
    distances[i] = dist_array[i];
  }

  return distances;
}

// gets the geometric mean of a vector of distances
//
// @param sites vector of distances
// @return geometric mean
// simulates de novos weighted by mutation rate
double _geomean(double distances[], int & len, double p, double offset) {
    /**
        gets the geometric mean of a vector of distances
        
        @sites vector of distances
        @return geometric mean
    */
    double total = 0;
    if(p != 0) {
      for (int i = 0; i< len; i++){
	total += std::pow(distances[i],p)/len;
      }
      return std::pow(total,1.0/p);
    } else {
      // if some values are zero, adjust the values upwards, then add the log10
      // value, otherwise add the uncorrected log10 value
      bool zero_val = _has_zero(distances, len);
      if (zero_val) {
        for (int i=0; i < len; i++) {
	  //total += log10(distances[i] + 1);
	  total += log10(distances[i] + offset);
        }
      } else {
        for (int i=0; i < len; i++) {
	  total += log10(distances[i]);
        }
      }
    
      // calculate the mean value
      double mean = total/len;
      mean = std::pow(10, mean);
    
      // adjust mean back to where it should be if we had a zero value
      if (zero_val) {
	//mean -= 1;
	mean -= offset;
      }
    
      return mean;
    }
}

double _geomean(std::vector<double> distances, double p, double offset) {
    double * dist = distances.data();
    int len = distances.size();
    return _geomean(dist, len, p, offset);
}

// double _geomean(double distances[], int & len, double p) {
//     /**
//         gets the geometric mean of a vector of distances
        
//         @sites vector of distances
//         @return geometric mean
//     */
//     double total = 0;
//     if( p != 0) {
//       for (int i = 0; i< len; i++){
// 	total += std::pow(distances[i],p)/len;
//       }
//       return std::pow(total,1.0/p);
//     } else {
//       // if some values are zero, adjust the values upwards, then add the log10
//       // value, otherwise add the uncorrected log10 value
//       bool zero_val = _has_zero(distances, len);
//       if (zero_val) {
//         for (int i=0; i < len; i++) {
// 	  total += log10(distances[i] + 1);
//         }
//       } else {
//         for (int i=0; i < len; i++) {
// 	  total += log10(distances[i]);
//         }
//       }
      
//       // calculate the mean value
//       double mean = total/len;
//       mean = std::pow(10, mean);
      
//       // adjust mean back to where it should be if we had a zero value
//       if (zero_val) { mean -= 1; }

//       return mean;
//     }
// }


// double _geomean(std::vector<double> distances, double p) {
//     double * dist = distances.data();
//     int len = distances.size();
//     return _geomean(dist, len, p);
// }

// @param choices Chooser object, to sample sites
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @return a list of mean distances for each iteration

std::vector<double> _simulate_distribution(Chooser & choices,
					   int iterations,
					   int de_novo_count,
					   double p,
					   std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    double distances[distance_len];
    int positions[de_novo_count];
    double scores[de_novo_count];
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
        // randomly select de novo sites for the iteration
      choices.generate_choices(positions, scores, de_novo_count);
        _get_distances(positions, de_novo_count, distances);

	if(scores[0] != -1)
	  _scale_distances(distances, scores, de_novo_count);
	mean_distances[n] = _geomean(distances, distance_len, p, 1);
    }
    
    // make sure the mean distances are sorted, so we can quickly merge with
    // previous distances
    std::sort(mean_distances.begin(), mean_distances.end());
    
    if(dist_file_output != ""){
      std::cout << "HERE IT IS " << dist_file_output << "\n" << std::endl;
      std::ofstream ofs(dist_file_output);
      std::ostream_iterator<double> output_iterator(ofs, "\n");
      std::copy(mean_distances.begin(), mean_distances.end(), output_iterator);
    }
    
    return mean_distances;
}


std::vector<double> _simulate_distribution_entropy(Chooser & choices,
						   int iterations,
						   int de_novo_count,
						   double p,
						   std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> entropies(iterations, 0);
    //int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    //double distances[distance_len];
    int positions[de_novo_count];
    double scores[de_novo_count];
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
      entropies[n] = 0;
      std::map<double,int> counts;
      // randomly select de novo sites for the iteration
      choices.generate_choices(positions, scores, de_novo_count);

      // count number of recurrences of each position
      for (int i=0; i < de_novo_count; i++) {
	positions[i] = choices.choice().pos;
	counts[floor((float)positions[i] / 3) + 1] += 1;
      }

      // calculate the entropy
      for (const auto &myPair : counts) {
	//LATER: offset by S0 entropy
	//entropies[n] += log(de_novo_count) - -1*myPair.second/(float)de_novo_count * log(myPair.second/(float)de_novo_count);
	entropies[n] += -1*myPair.second/(float)de_novo_count * log(myPair.second/(float)de_novo_count);	
      }
    }
    
    // make sure the entropies are sorted, so we can quickly merge with
    // previous distances
    std::sort(entropies.begin(), entropies.end());
    
    return entropies;
}


std::vector<double> _simulate_distribution(Chooser & choices,
					   double three_d_locations[][3],
					   int iterations,
					   int de_novo_count,
					   double p,
					   std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */

    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    int positions[de_novo_count];
    double distances[distance_len];
    double scores[de_novo_count];
    // run through the required iterations

    std::cout << "starting with" << iterations <<" iterations" << std::endl;
    for (int n=0; n < iterations; n++) {
      // randomly select de novo sites for the iteration
      choices.generate_choices(positions, scores, de_novo_count);
      _get_distances(positions, three_d_locations, de_novo_count, distances);
      if(scores[0] != -1)
	_scale_distances(distances, scores, de_novo_count);
      mean_distances[n] = _geomean(distances, distance_len, p, 3.5/*1*/);
    }

    std::sort(mean_distances.begin(), mean_distances.end());
    if(dist_file_output != ""){
      std::cout << "HERE IT IS" << dist_file_output << "\n" << std::endl;
      std::ofstream ofs(dist_file_output);
      std::ostream_iterator<double> output_iterator(ofs, "\n");
      std::copy(mean_distances.begin(), mean_distances.end(), output_iterator);
    }

    
    return mean_distances;
}



std::vector<double> _simulate_distribution_coevol(Chooser & choices,
						  const std::vector<double> coevol,
						  int iterations,
						  int de_novo_count,
						  double p,
						  std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    double distances[distance_len];
    int positions[de_novo_count];
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
      // randomly select de novo sites for the iteration
      choices.generate_choices(positions, de_novo_count);
      // convert the positions into distances between all pairs, and get the
      // geometric mean distance of all the distances
      _get_distances_coevol(positions,
			    coevol.data(),
			    sqrt(coevol.size()),
			    de_novo_count,
			    distances);
        
      mean_distances[n] = _geomean(distances, distance_len, p, 1);
    }
    // make sure the mean distances are sorted, so we can quickly merge with
    // previous distances
    std::sort(mean_distances.begin(), mean_distances.end());
    if(dist_file_output != ""){
      std::ofstream ofs(dist_file_output);
      std::ostream_iterator<double> output_iterator(ofs, "\n");
      std::copy(mean_distances.begin(), mean_distances.end(), output_iterator);
    }

    return mean_distances;
}


// halt permutations if the P value could never be significant.
//
// assess whether the P-value could never fall below 0.1, and cut out after a
// smaller number of iterations, in order to minimise run time. Figure out the
// lower bound of the confidence interval for the current simulated P value.
//
// @param p_val current simulated P value
// @param iterations iterations run in order to obtain the simulated P value
// @param z standard normal deviate (eg 1.96 for 95% CI)
// @param alpha value above which to cease permuting
// @return True/False for whether to halt the permuations
bool _halt_permutation(double p_val, int iterations, double z, double alpha) {
    // @todo figure out whether this is legit, perhaps a Sequential
    // @todo Probability Ratio Test (SPRT) would be more appropriate.

    double delta = (z * sqrt((p_val * (1 - p_val))))/iterations;
    double lower_bound = p_val - delta;
    
    // if the lower bound of the confidence interval exceeds 0.1, then we
    // can be sure it's not going to ever get lower than 0.05.
    return lower_bound > alpha;
}

// simulates de novos weighted by mutation rate
//
// @param choices Chooser object, to sample sites
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @param observed_value mean distance observed in the real de novo events
// @return a list of mean distances for each iteration
double _analyse_de_novos(Chooser & choices,
			 int iterations,
			 int de_novo_count,
			 double observed_value,
			 double p,
			 std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;
    
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution(choices,
							      iters_to_run,
							      de_novo_count,
							      p,
							      dist_file_output);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());

	// halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        
        iterations += 1000000;  // for if we need to run more iterations
    }
    
    return sim_prob;
}


double _analyse_de_novos_entropy(Chooser & choices,
				 int iterations,
				 int de_novo_count,
				 double observed_value,
				 double p,
				 std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;
    
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution_entropy(choices,
							      iters_to_run,
							      de_novo_count,
							      p,
							      dist_file_output);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());

	// halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        
        iterations += 1000000;  // for if we need to run more iterations
    }
    
    return sim_prob;
}

double _analyse_de_novos(Chooser & choices,
			 std::vector<std::vector<double>> three_d_locations,
			 int iterations,
			 int de_novo_count,
			 double observed_value,
			 double p,
			 std::string dist_file_output)
{
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;

    int aa_count = three_d_locations.size();
    double three_d_locations_array[aa_count][3];
    for (int i=0; i<aa_count; i++) {
        three_d_locations_array[i][0] = three_d_locations[i][0];
        three_d_locations_array[i][1] = three_d_locations[i][1];
        three_d_locations_array[i][2] = three_d_locations[i][2];	
    }
    std::cout << "Starting simulation" << std::endl;
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));

        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution(choices,
							      three_d_locations_array,
							      iters_to_run,
							      de_novo_count,
							      p,
							      dist_file_output);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());
	
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        iterations = iterations*2;
        //iterations += 1000000;  // for if we need to run more iterations
    }

    std::cout << "Sim prob = " << sim_prob << std::endl;
    return sim_prob;
}

std::vector<double> _simulate_distribution_multi(std::vector<std::string> chains,
						 std::vector<std::string> proteins,
						 std::vector<int> de_novo_count_per_chain,
						 std::map<std::string,Chooser> choices,						 
						 //std::vector<Chooser> choices,
						 //						 std::vector<string> proteins_of_choices,
						 std::vector<std::map<int,std::vector<double>>> three_d_locations,
						 int iterations,
						 int de_novo_count,
						 double p,
						 std::string dist_file_output) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
  
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    //    int positions[de_novo_count];
    double distances[distance_len];
    //    double scores[de_novo_count];
    std::vector<double> scores;
    // run through the required iterations
    int m = 0;
    for (int n=0; n < iterations; n++) {
      // randomly select de novo sites for the iteration
      std::map<int, std::vector<int>> positions;
      std::map<int, std::vector<double>> scores_vector;      
      std::vector<double> xyz_positions;
      std::vector<int> positions_rec;
      m = 0;
      std::vector<std::string> seen_proteins;
      for( auto c : chains) {
	//	std::cout << "chain is " << c << " and protein is " << proteins[m] << " and it has " << de_novo_count_per_chain[m] << " de novos"<< std::endl;

	if( de_novo_count_per_chain[m] == 0 )
	  continue;
	
	// Is this current chain of a protein we have already seen?
	// If so, do not regenerate random choices but use the
	// same ones generated before.
	bool protein_not_seen_before = (m == 0) ? true : std::count(seen_proteins.begin(), seen_proteins.end(), proteins[m])==0;

	if(protein_not_seen_before) {
	  //generate the random choices
	  choices[proteins[m]].generate_choices(positions[m],
						scores_vector[m],
						de_novo_count_per_chain[m]);
	  //for(auto p : positions[m])
	  //  std::cout << p << " ";
	  //std::cout <<std::endl;
	  seen_proteins.push_back(proteins[m]);
	} else {
	  //recycle choices on the chain
	  auto first_occurence_idx = std::find(chains.begin(), std::next(chains.begin(),m + 1), chains[m]);
	  int first_occurence = std::distance(chains.begin(), first_occurence_idx)-1;

	  positions[m] = positions[first_occurence];
	  scores_vector[m] = scores_vector[first_occurence];
	  //assert(chains[m] == chains[first_occurence]);
	  //assert(positions[m].size() == de_novo_count_per_chains[m]);
	}

	//Set three d locations for current chain and add to list
	int codon;
	for (auto p : positions[m]){
	  codon = floor((float)p / 3) + 1;
	  positions_rec.push_back(codon);
	  for( auto i : three_d_locations[m][codon]){
	    xyz_positions.push_back(i);
	  }
	}
	
	//Add scores at proper index
	for (auto score : scores_vector[m]){
	  scores.push_back(score);
	}
	m = m + 1;
      }
      
      //assert(xyz_positions.size() == de_novo_count*3);
      //get distances between all de novos
      int idx = 0;
      for(int i = 0; i < de_novo_count; i++){
	for(int j = i+1; j < de_novo_count; j++){
	  distances[idx] = 
	    std::pow((xyz_positions[3*i]-xyz_positions[3*j])*(xyz_positions[3*i]-xyz_positions[3*j]) +
		   (xyz_positions[3*i+1]-xyz_positions[3*j+1])*(xyz_positions[3*i+1]-xyz_positions[3*j+1]) +
		   (xyz_positions[3*i+2]-xyz_positions[3*j+2])*(xyz_positions[3*i+2]-xyz_positions[3*j+2]), .5);
	  idx = idx + 1;
	}
      }
      //assert(idx == distance_len);
      if(scores[0] != -1)
	_scale_distances(distances, scores.data(), de_novo_count);
      mean_distances[n] = _geomean(distances, distance_len, p, 3.5/*1 multi*/);

      if(std::isnan(mean_distances[n]) || std::isinf(mean_distances[n])){
	std::cout << "HERE" << std::endl;
	for(int it = 0; it < distance_len; it++){
	  std::cout << distances[it] << " ";
	}

	for(int it = 0; it < de_novo_count; it++){
	  std::cout << xyz_positions[3*it] << " " <<xyz_positions[3*it + 1] << " " <<xyz_positions[3*it+2]<< " " 
		    << positions_rec[it] << std::endl;
	}
	
      }
      
    }

    std::sort(mean_distances.begin(), mean_distances.end());


    int l = mean_distances.size();
    std::cout <<l << std::endl;
    std::cout << mean_distances[floor(l/10)]<< std::endl;
    std::cout << mean_distances[floor(l*2/10)]<< std::endl;
    std::cout << mean_distances[floor(l*3/10)]<< std::endl;
    std::cout << mean_distances[floor(l*4/10)]<< std::endl;
    std::cout << mean_distances[floor(l*5/10)]<< std::endl;
    std::cout << mean_distances[floor(l*6/10)]<< std::endl;
    std::cout << mean_distances[floor(l*7/10)]<< std::endl;
    std::cout << mean_distances[floor(l*8/10)]<< std::endl;
    std::cout << mean_distances[floor(l*9/10)]<< std::endl;       std::cout << mean_distances[floor(l)-1]<< std::endl;

    if(dist_file_output != ""){
      std::ofstream ofs(dist_file_output);
      std::ostream_iterator<double> output_iterator(ofs, "\n");
      std::copy(mean_distances.begin(), mean_distances.end(), output_iterator);
    }

    return mean_distances;
}


double _analyse_de_novos_multi(std::vector<std::string> chains,
			       std::vector<std::string> proteins_of_chains,
			       std::vector<int> de_novo_count_per_chain,			       
			       std::map<std::string,Chooser> choices,
			       //std::vector<Chooser> choices,
			       //			       std::vector<std::string> proteins_of_choices,
			       std::vector<std::map<int,std::vector<double>>> three_d_locations,
			       int iterations,
			       int de_novo_count,
			       double observed_value,
			       double p,
			       std::string dist_file_output)
{
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;

    //for(auto c : three_d_locations){
    //  for(auto v : c) {
    //	std::cout << v.first << std::endl;
    //  std::cout << v.second[0] << " " << v.second[1] << " " << v.second[2] << std::endl;
    //  }
    // }

    
    // int aa_count = three_d_locations.size();
    //double three_d_locations_array[aa_count][3];
    //for (int i=0; i<aa_count; i++) {
    //    three_d_locations_array[i][0] = three_d_locations[i][0];
    //    three_d_locations_array[i][1] = three_d_locations[i][1];
    //    three_d_locations_array[i][2] = three_d_locations[i][2];	
    //}
    std::cout << "Starting simulation" << std::endl;
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));

        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution_multi(chains,
								    proteins_of_chains,
								    de_novo_count_per_chain,
								    choices,
								    //proteins_of_choices,
								    three_d_locations,
								    iters_to_run,
								    de_novo_count,
								    p,
								    dist_file_output);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());
	
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        iterations = iterations*2;
        //iterations += 1000000;  // for if we need to run more iterations
    }

    std::cout << "Sim prob = " << sim_prob << std::endl;
    return sim_prob;
}



double _analyse_de_novos_coevol(Chooser & choices,
				std::vector<double> coevol,
				int iterations,
				int de_novo_count,
				double observed_value,
				double p,
				std::string dist_file_output)
{
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;

    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution_coevol(choices,
								     coevol,
								     iters_to_run,
								     de_novo_count,
								     p,
								     dist_file_output);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());
	sim_prob = 1.0 - sim_prob;
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        
        iterations += 1000000;  // for if we need to run more iterations
    }

    std::cout << "Sim prob = " << sim_prob << std::endl;
    return sim_prob;
}

double get_distance(int i,
		    int j,
		    int positions [],
		    std::vector<std::vector<double>> xyz){

  int codon_i = floor((float)positions[i] / 3);
  int codon_j = floor((float)positions[j] / 3);  

  return std::pow((xyz[codon_i][0] - xyz[codon_j][0])*(xyz[codon_i][0] - xyz[codon_j][0]) +
		  (xyz[codon_i][1] - xyz[codon_j][1])*(xyz[codon_i][1] - xyz[codon_j][1]) +
		  (xyz[codon_i][2] - xyz[codon_j][2])*(xyz[codon_i][2] - xyz[codon_j][2]), .5);
}

std::vector<double> _simulate_distribution_west(Chooser & choices,
						std::vector<std::vector<double>> three_d_locations,
						int iterations,
						int de_novo_count) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */

    // use a vector to return the mean distances, easier to call from python
    std::vector<double> compiled_scores(iterations);
    int positions[de_novo_count];
    double scores[de_novo_count];
    double total_score = 0;
    // run through the required iterations

    std::cout << "starting with" << iterations <<" iterations" << std::endl;

    bool scaled = false;
    for (int n=0; n < iterations; n++) {
      // randomly select de novo sites for the iteration
      choices.generate_choices(positions, scores, de_novo_count);
      if (scaled == false){
	total_score = std::accumulate(scores,
				      scores+de_novo_count,
				      0.0);
      }
      else{
	double new_scores[de_novo_count];
	double distance = 0;
	for(int i = 0; i < de_novo_count; i++){
	  for(int j = 0; j < de_novo_count; j++){
	    if( i < j){
	      distance = get_distance(i,j,positions, three_d_locations);
	      new_scores[i] += scores[j]*distance;
	      new_scores[j] += scores[j]*distance;
	    }	   
	    total_score += scores[i];
	  }
	}
      }
      compiled_scores.push_back(total_score);
    }
    std::sort(compiled_scores.begin(), compiled_scores.end());    
    return compiled_scores;
}
    
double prob_score_with_k_denovos(double observed_value,
				 int de_novo_count,
				 Chooser & choices,
				 std::vector<std::vector<double>> three_d_locations,
				 int iterations)
{
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;

    std::cout << "Starting simulation" << std::endl;
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution_west(choices,
								   three_d_locations,
								   iters_to_run,
								   de_novo_count);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());
	
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        iterations = iterations*2;
        //iterations += 1000000;  // for if we need to run more iterations
    }

    std::cout << "Sim prob = " << sim_prob << std::endl;
    return sim_prob;
}

double prob_of_k_denovos(int k, double lambda){

  return pow(M_E, k * log(lambda) - lambda - lgamma(k + 1.0));
  
}

double _analyse_de_novos_west(Chooser & choices,
			      std::vector<std::vector<double>> three_d_locations,
			      int iterations,
			      double lambda,
			      int de_novo_count,
			      double observed_value,
			      double p)
{
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */

    double prob = 0;
    //double current_prob;
    for(int k  = 0; k<250; k++){
      prob += prob_of_k_denovos(k, lambda) * prob_score_with_k_denovos(observed_value,
								       k,
								       choices,
								       three_d_locations,
								       iterations);
    }

    return prob;
}
