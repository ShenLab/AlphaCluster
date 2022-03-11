#ifndef DENOVONEAR_SITESCHECKS_H_
#define DENOVONEAR_SITESCHECKS_H_

#include <string>
#include <vector>
#include <map>

#include "tx.h"
#include "weighted_choice.h"

class SitesChecks {
    /**
    class to build weighted choice random samplers for nonsense, missense,
    and functional classes of variants, using site specific mutation rates
    
    Only include the site specific probability if it mutates to a different
    amino acid, or occurs close to an intron/exon boundary, using consequences
    defined at: http://www.ensembl.org/info/genome/variation/predicted_data.html
    */
    
    std::unordered_map<std::string, std::unordered_map<std::string, double>> mut_dict;
    std::unordered_map<std::string, Chooser> rates;
    std::map<int, std::map<char, double>> scores;
    int boundary_dist;
    int kmer_length;
    int mid_pos;
    
    std::unordered_map<char, char> transdict = {
        {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};
    
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    std::vector<std::string> categories = {"missense", "nonsense", "synonymous",
        "splice_lof", "splice_region", "loss_of_function", "intronic"};

 public:
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords) :
         _tx { tx }, use_cds_coords { cds_coords } { init(mut); };
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords, int start, int end) :
        _tx { tx }, use_cds_coords { cds_coords } { init(mut, start, end); };
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords, std::vector<int> variants) :
     _tx { tx }, use_cds_coords { cds_coords } { init(mut, variants);};
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords,
	      std::map<int,std::map<char,double>> scores) :
    _tx { tx }, use_cds_coords { cds_coords } { init(mut, scores);};
    SitesChecks(Tx tx,
		std::vector<std::vector<std::string>> mut,
		bool cds_coords,
		std::map<int,std::map<char,double>> scores,
		std::vector<int> residues) :
      _tx { tx }, use_cds_coords { cds_coords } { init(mut, scores,residues);};  
    SitesChecks(Tx tx,
		std::vector<std::vector<std::string>> mut,
		bool cds_coords,
		std::map<int,std::map<char,double>> scores,
		double threshold) :
      _tx { tx }, use_cds_coords { cds_coords } { init(mut, scores, threshold);};  

    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords, Tx mask) :
         _tx { tx }, masked { mask }, use_cds_coords { cds_coords } { has_mask = true; init(mut); };
    Chooser * __getitem__(std::string category) { return &rates[category]; };
    void initialise_choices();
    
    Tx _tx;
    void set_inherited_controls(std::vector<int> inherited_variants, std::string category);
    void check_position(int bp, double threshold = -2);
    double check_score(char initial_aa, char mutated_aa, int position);
    void print_all(std::string category);
    int get_offset(int bp);
    void check_consequence(std::string & cq, char & initial_aa, char & mutated_aa, int & offset);
    
 private:
    Tx masked = Tx("zz", "z", -100, -100, '+');
    void init(std::vector<std::vector<std::string>> mut);
    void init(std::vector<std::vector<std::string>> mut, int start, int end);
    void init(std::vector<std::vector<std::string>> mut, std::vector<int> variants);
    void init(std::vector<std::vector<std::string>> mut, std::map<int,std::map<char,double>> scores);
    void init(std::vector<std::vector<std::string>> mut,
	      std::map<int,std::map<char,double>> scores,
	      std::vector<int> residues);
    void init(std::vector<std::vector<std::string>> mut,
	      std::map<int,std::map<char,double>> scores,
	      double threshold);  
    bool has_mask = false;
    bool use_cds_coords = true;
    bool scores_set = false;
};
Region _get_gene_range(Tx & tx);
char _get_mutated_aa(Tx & tx, char base, std::string codon, int intra_codon);

#endif  // DENOVONEAR_SITESCHECKS_H_
