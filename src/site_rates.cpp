#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include <iostream>

#include "site_rates.h"

// get the lowest and highest positions of a transcripts coding sequence
Region _get_gene_range(Tx & tx) {
    int boundary_1 = tx.get_cds_start();
    int boundary_2 = tx.get_cds_end();
    
    int start = std::min(boundary_1, boundary_2);
    int end = std::max(boundary_1, boundary_2);
    
    return Region {start, end};
}

// find the amino acid resulting from a base change to a codon
//
// @param tx transcript object for a gene
// @param base alternate base (e.g. 'G') to introduce
// @param codon DNA sequence of a single codon
// @param intra_codon position within the codon to be altered (0-based)
// @returns single character amino acid code translated from the altered
//     codon.
std::string _get_mutated_aa(Tx & tx, char base, std::string codon, int intra_codon) {
    // figure out what the mutated codon is
    codon[intra_codon] = base;
    return tx.translate(codon);
}

void SitesChecks::init(std::vector<std::vector<std::string>> mut) {
    
    // convert the rates data to a nested map
    for (auto line : mut) {
        mut_dict[line[0]][line[1]] = std::stod(line[2]);
    }
    
    // calculate the length of sequences used in the mutation rate dictionary
    // This means we can flexibly use 3-mers, 5-mers, 7-mers etc if desired.
    std::string str = mut_dict.begin()->first;
    kmer_length = str.length();
    mid_pos = kmer_length/2;
    
    initialise_choices();
    // check the consequence alternates for each base in the coding sequence
    Region region = _get_gene_range(_tx);
    for (int i=region.start; i < region.end + 1; i++ ) {
        check_position(i);
    }
}

// initialise a WeightedChoice object for each consequence category
void SitesChecks::initialise_choices() {
    for (auto category : categories) {
        rates[category] = Chooser();
    }
}

// get the consequence of an amino acid change (or not)
std::string SitesChecks::check_consequence(std::string & initial_aa,
        std::string & mutated_aa, int & position) {
    std::string cq = "synonymous";
    
    bool in_coding = _tx.in_coding_region(position);
    if ( initial_aa != "*" && mutated_aa == "*" ) {
        // checks if two amino acids are a nonsense (eg stop_gained) mutation
        cq = "nonsense";
    } else if ( !in_coding && boundary_dist < 3 ) {
        // check if a variant has a splice_donor or splice_acceptor consequence
        // These variants are defined as being the two intronic base-pairs
        // adjacent to the intron/exon boundary.
        cq = "splice_lof";
    } else if ( initial_aa != mutated_aa ) {
        // include the site if it mutates to a different amino acid.
        cq = "missense";
    } else if (in_coding) {
        if ( boundary_dist < 4) {
          // catch splice region variants within the exon, and in the appropriate
          // region of the intron (note that loss of function splice_donor and
          // splice_acceptor variants have been excluded when we spotted nonsense).
          cq = "splice_region";
        }
    } else {
        if (boundary_dist < 9) {
            // check for splice_region_variant inside intron
            cq = "splice_region";
        } else {
            cq = "intronic";
        }
    }
    
    return cq;
}

// add the consequence specific rates for the alternates for a variant
//
// @param bp genomic position of the variant
void SitesChecks::check_position(int bp) {
    if (bp == 2684742) { std::cout << "checking position: " << bp << "\n"; }
    // ignore sites within masked regions (typically masked because the
    // site has been picked up on alternative transcript)
    if ( has_mask && masked.in_coding_region(bp) ) {
        return ;
    }
    
    // ignore sites outside the CDS region
    if (bp < std::min(_tx.get_cds_start(), _tx.get_cds_end()) ||
        bp > std::max(_tx.get_cds_start(), _tx.get_cds_end())) {
        return ;
    }
    
    if (bp == 2684742) { std::cout << " - getting seq and distance\n";}
    std::string seq = _tx.get_centered_sequence(bp, kmer_length);
    boundary_dist = _tx.get_boundary_distance(bp);
    
    char fwd = '+';
    if (_tx.get_strand() != fwd) {
        seq = _tx.reverse_complement(seq);
    }
    
    if (bp == 2684742) { std::cout << " - getting codon info\n";}
    Codon codon;
    try {
        codon = _tx.get_codon_info(bp);
    } catch ( const std::invalid_argument& e ) {
        return ;
    }
    
    std::string initial_aa = codon.initial_aa;
    int cds_pos = codon.cds_pos;
    int offset = codon.offset;

    if (bp == 2684742) { std::cout << " - assigning alts\n";}
    // std::cout << " - assigning alts\n";
    // drop the initial base, since we want to mutate to other bases
    std::vector<char> alts(bases);
    if (bp == 2684742) { std::cout << " - dropping ref from alts\n";}
    char ref = seq[mid_pos];
    if (bp == 2684742) { std::cout << "initial_aa: " << initial_aa << "\n";}
    if (bp == 2684742) { std::cout << "boundary_dist: " << boundary_dist << "\n"; }
    if (bp == 2684742) { std::cout << "seq: " << seq << "\n";}
    if (bp == 2684742) { std::cout << "ref: " << ref << "\n";}
    if (bp == 2684742) { std::cout << "alts: " << alts[0] << alts[1] << alts[2] << alts[3] << "\n";}
    alts.erase(std::find(alts.begin(), alts.end(), ref));
    
    if (bp == 2684742) { std::cout << " - preparing variables\n";}
    std::string mutated_aa;
    std::string alt_seq = seq;
    std::string category;
    for (auto &alt : alts) {
        mutated_aa = initial_aa;
        // alt_seq = seq.substr(0, mid_pos) + alt +
        //     seq.substr(mid_pos + 1, kmer_length - mid_pos);
        if (bp == 2684742) {std::cout << " - - checking alt: " << alt <<"\n";}
        
        alt_seq[mid_pos] = alt;
        
        // std::cout << " - - getting mutated AA\n";
        double rate = mut_dict[seq][alt_seq];
        if ( initial_aa.size() != 0 ) {
            mutated_aa = _get_mutated_aa(_tx, alt, codon.codon_seq, codon.intra_codon);
        }
        
        // std::cout << " - - checking consequence\n";
        category = check_consequence(initial_aa, mutated_aa, bp);
        
        // figure out what the ref and alt alleles are, with respect to
        // the + strand.
        // std::cout << " - - flipping alleles\n";
        if (_tx.get_strand() != fwd) {
            ref = transdict[ref];
            alt = transdict[alt];
        }
        
        // std::cout << " - - adding choice\n";
        if (use_cds_coords) {
            rates[category].add_choice(cds_pos, rate, ref, alt, offset);
        } else {
            rates[category].add_choice(bp, rate, ref, alt, 0);
        }
        
        // std::cout << " - - checking lof\n";
        if (category == "nonsense" || category == "splice_lof") {
            rates["loss_of_function"].add_choice(cds_pos, rate, ref, alt, offset);
        }
    }
}
