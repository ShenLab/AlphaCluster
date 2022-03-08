#!/usr/bin/env python

import torch
from torch import nn 
import pickle
import argparse
import numpy as np

def run(args):
    # Pick transcript
    transcript_id = args.transcript

    # Fetch file
    feature_dir = args.feature_dir
    transcript_feature_file = f'{feature_dir}/{transcript_id}.pickle'
    with open(transcript_feature_file, 'rb') as infile:
        transcript_features = pickle.load(infile)
        feature_len = transcript_features.shape[0]

    # 20 AA and gap
    A = 21
    # How many species to use for the calculation
    # Maximum is 200
    num_species = 200

    # Grab all available species and limit to num_species
    x = transcript_features[:,21:421]
    x = x[:,:num_species]
    x = torch.tensor(x)
        
    # One hot encode, dim is sites_i * species_s * AA_a
    x = nn.functional.one_hot(x.type(torch.int64), num_classes=21 ).type(torch.float32)

    # Create one site Kronecker delta function
    # dim is sites_i x AA_a
    # entry [ia] is freq of occurences of AA_a at site_i across species
    one_site_freq = 1/num_species * torch.einsum('isa->ia', x)
        
    # Create a dual site Kronecker delta function
    # dim is sites_i * sites_j * AA_a * AA_b
    # entry [ijab] is  freq of occurences of AA_a at site_i * freq of occurrences of AA_b at site_j
    dual_one_site_freq = torch.einsum('ia,jb->ijab',
                                          one_site_freq, one_site_freq)
        
    # Create two site Kronecker delta function 
    # dim is sites_i * sites_j * AA_a * AA_b
    # entry [ijab] is sum of (freq of simulanteious occurences of AA_a at site_i and AA_b at site_j)
    two_site_freq = 1/num_species * torch.einsum('isa,jsb->ijab',x,x)
        
    # Calculate covariance
    # A 21 x 21 matrix for each combination of site_i and site_j
    coev_corr_per_aa_per_site = (two_site_freq - dual_one_site_freq)

    # Calculate the correlation for a given site pair ij
    coev_corr_per_site = torch.sqrt(
        torch.einsum('ijab,ijab->ij',
                     coev_corr_per_aa_per_site,
                     coev_corr_per_aa_per_site)
    )
        
    np.savetxt(args.output, coev_corr_per_site.numpy(),fmt='%1.5f') 


def main():
    msg = "Calculate coevolution strength for a gene from MSAs"

    # Initialize parser
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("-t", "--transcript", help = "Name transcript. For example, 'ENST00000325103'", required=True)
    parser.add_argument("-f", "--feature_dir", help = "gMVP feature directory", default = "/share/terra/backups/MVPContext/combined_feature_2021_v2")
    parser.add_argument("-o", "--output", help = "Name of output file", required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main()
