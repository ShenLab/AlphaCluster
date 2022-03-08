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
void _get_distances_coevol(const std::vector<int> & sites,
			   const std::vector<double> & coevol,
			   std::vector<double> & distances);
void _get_distances_coevol(const int sites[],
			   const double coevol[],
			   const int aa_length,
			   int & len,
			   double distances[]);

void _get_distances(int sites[], 
		    int & len, 
		    double distances[]);
void _get_distances(const std::vector<int> & sites, 
		    std::vector<double> & distances);

void _scale_distances(double distances[],
		      const double scores[],
		      int len);

std::vector<double> _scale_distances(std::vector<double> distances,
		      const std::vector<double> scores);

bool _has_zero(double distances[], 
	       int & len);
bool _has_zero(std::vector<double> distances);
//bool _has_zero(double distances[], 
//	       int & len);
//bool _has_zero(std::vector<double> distances);
double _geomean(double distances[], 
		int & len, double p, double offset);
double _geomean(std::vector<double> distances, double p, double offset);
//double _geomean(int distances[], 
//		int & len, double p);
//double _geomean(std::vector<int> distances, double p);
bool _halt_permutation(double p_val, 
		       int iterations, 
		       double z = 10.0,
		       double alpha = 0.01);
std::vector<double> _simulate_distribution(Chooser & choices,
					   int iterations, 
					   int de_novo_count,
					   double p,
					   std::string dist_file_output);
std::vector<double> _simulate_distribution_entropy(Chooser & choices,
						   int iterations, 
						   int de_novo_count,
						   double p,
						   std::string dist_file_output);

std::vector<double> _simulate_distribution(Chooser & choices,
					   double three_d_locations[][3],
					   int iterations, 
					   int de_novo_count,
					   double p,
					   std::string dist_file_output);
std::vector<double> _simulate_distribution_coevol(Chooser & choices,
						  std::vector<double> & coevol,
						  int iterations, 
						  int de_novo_count,
						  double p,
						  std::string dist_file_output);

double _analyse_de_novos(Chooser & choices, 
			 int iterations,
			 int de_novo_count,
			 double observed_value,
			 double p,
			 std::string dist_file_output);

double _analyse_de_novos_entropy(Chooser & choices, 
				 int iterations,
				 int de_novo_count,
				 double observed_value,
				 double p,
				 std::string dist_file_output);

double _analyse_de_novos(Chooser & choices,
			 std::vector<std::vector<double>> three_d_locations, 
			 int iterations,
			 int de_novo_count,
			 double observed_value,
			 double p,
			 std::string dist_file_output);
double _analyse_de_novos_coevol(Chooser & choices,
				std::vector<double> coevol, 
				int iterations,
				int de_novo_count,
				double observed_value,
				double p,
				std::string dist_file_output);

double _analyse_de_novos_west(Chooser & choices,
			      std::vector<std::vector<double>> three_d_locations,
			      int iterations,
			      double lambda,
			      int de_novo_count,
			      double observed_value,
			      double p);

double _analyse_de_novos_multi(std::vector<std::string> chains,
			       std::vector<std::string> proteins,			       
			       std::vector<int> de_novo_count_per_chain,
			       std::map<std::string,Chooser> choices,
			       //std::vector<Chooser> choices,
			       //std::vector<std:string> proteins_of_choices,
			       std::vector<std::map<int,std::vector<double>>> three_d_locations,
			       int iterations,
			       int de_novo_count,
			       double observed_value,
			       double p,
			       std::string dist_file_output);
std::vector<double> _simulate_distribution_multi(std::vector<std::string> chains,
						 std::vector<std::string> proteins_of_chains,
						 std::vector<int> de_novo_count_per_chain,
						 std::map<std::string,Chooser> choices,						 
						 //std::vector<Chooser> choices,
						 //std::vector<std:string> proteins_of_choices,
						 std::vector<std::map<int,std::vector<double>>> three_d_locations,
						 int iterations,
						 int de_novo_count,
						 double p,
						 std::string dist_file_output);
#endif  // DENOVONEAR_SIMULATE_H_
