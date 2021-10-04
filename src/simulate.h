#ifndef DENOVONEAR_SIMULATE_H_
#define DENOVONEAR_SIMULATE_H_

#include <vector>
#include <map>


void _get_distances(int sites[],
		    double three_d_locations[][3],
		    //const std::vector<std::vector<double>> & three_d_locations,
		    int & len,
		    double distances[]);
void _get_distances(const std::vector<int> & sites,
		    //double three_d_locations[][3],		    
		    const std::vector<std::vector<double>> & three_d_locations,
		    std::vector<double> & distances);
void _get_distances(int sites[], 
		    int & len, 
		    int distances[]);
void _get_distances(const std::vector<int> & sites, 
		    std::vector<int> & distances);
bool _has_zero(int distances[], 
	       int & len);
bool _has_zero(std::vector<int> distances);
bool _has_zero(double distances[], 
	       int & len);
bool _has_zero(std::vector<double> distances);
double _geomean(double distances[], 
		int & len);
double _geomean(std::vector<double> distances);
double _geomean(int distances[], 
		int & len);
double _geomean(std::vector<int> distances);
bool _halt_permutation(double p_val, 
		       int iterations, 
		       double z = 10.0,
		       double alpha = 0.01);
std::vector<double> _simulate_distribution(Chooser & choices,
					   int iterations, 
					   int de_novo_count);
std::vector<double> _simulate_distribution(Chooser & choices,
					   double three_d_locations[][3],
					   //const std::vector<std::vector<double>> & three_d_locations,
					   int iterations, 
					   int de_novo_count);
double _analyse_de_novos(Chooser & choices, 
			 int iterations,
    int de_novo_count, double observed_value);
double _analyse_de_novos(Chooser & choices,
			 //double three_d_locations[][],			 
			 std::vector<std::vector<double>> three_d_locations, 
			 int iterations,
			 int de_novo_count, double observed_value);

#endif  // DENOVONEAR_SIMULATE_H_
