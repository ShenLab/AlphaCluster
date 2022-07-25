import subprocess
import os
import pipes
import pysam.bcftools
from pysam import VariantFile
import numpy as np

def load_inherited(chrom,
                   start,
                   end,
                   conseqs = ["missense_variant"],
                   max_af = 0.0001,
                   setting = "Affected"):

    # Import phenotype data from SPARK 
    affected = []
    unaffected = []

    mastertables = ["/share/terra/SPARK_WES_1/Mastertables/SPARK.WES1.mastertable.2021_03.tsv",
                    "/share/terra/SPARK_WES_2/Mastertables/SPARK.WES2.mastertable.2020_06.tsv"]

    for mastertable in mastertables:
        data = [x.split('\t') for x in open(mastertable).read().splitlines()]
        affected += [x[0] for x in data if x[5] == "2"]
        unaffected += [x[0] for x in data if x[5] == "1"]
    affected = set(affected)
    unaffected = set(unaffected)

    # Import variants from SPARK
    variants = []
    vcfs = ["/share/terra/SPARK_WES_1/Variants/DeepVariant/pvcf/wes1_27281_exome.ann.vcf.gz",
            "/share/terra/SPARK_WES_2/Variants/DeepVariant/pvcf/wes2_15995_exome.ann.vcf.gz"]

    for vcf_in in vcfs :
        for conseq in conseqs:
#        print("filtering on " + vcf_in + "for " + str(chrom)+":"+str(start)+"-"+str(end) )

            cmd_line = "bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%TGT\t%ALT\t%INFO/Consequence\n]' -r"  + str(chrom)+":"+str(start)+"-"+str(end) + " -i '(gnomAD_AF < " + str(max_af) + " || gnomAD_AF = " + '".")' + ' && GT = "alt" && INFO/Consequence == "' + str(conseq) +'" ' + "' " + vcf_in
        
            po = subprocess.Popen(cmd_line,
                                  shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

            stdout, stderr = po.communicate()

            po.wait()

            return_code = po.returncode

            if return_code != 0:
                raise Exception("Command line {} got return code {}.\nSTDOUT: {}\nSTDERR: {}".format(cmd_line, return_code, stdout, stderr))
            variants =  variants + [x.split('\t') for x in stdout.decode().split('\n')[:-1]]

    if setting == "Affected":
        variants = [(int(x[1]), x[5]) for x in variants if x[2] in affected]
    elif setting == "Unaffected":
        variants = [(int(x[1]),x[5]) for x in variants if x[2] in unaffected]

    return variants

def load_inherited_controls(chrom,
                            start,
                            end,
                            conseqs = ["missense_variant"],
                            max_af = 0.0001):
    
    # Locate approriate vcf
    pvcf_blocks_file = "/share/terra/UKBB/pVCF/pvcf_blocks.txt"
    start_block = 0
    end_block = 0
    with open(pvcf_blocks_file, 'r') as file:
        found_start = False
        for line in file:
            line = line.rstrip().split()
            if not int(line[1]) == int(chrom):
                continue
            if not found_start and int(line[3]) <= int(start) and int(line[4]) >= int(start):
                start_block = int(line[2])
                found_start = True
                
            if found_start and int(line[4]) >= int(end):
                end_block = int(line[2])
                break

    vcfs = []
    for vcf_block in range(start_block, end_block + 1):
        vcfs += ["/share/terra/UKBB/pVCF/ukb23156_c" + str(chrom) + "_b" + str(vcf_block) + "_v1.vcf.gz"]
        
    # Import variants from UKBB
    variants = []

    for vcf_in in vcfs :
        for conseq in conseqs:
            print("filtering on " + vcf_in + "for " + str(chrom) + ":" + str(start) + "-" + str(end))

            cmd_line = "bcftools view -M2 -m2 -v snps -r chr"+ str(chrom)+":"+str(start)+"-"+str(end) + " -i 'AF<=" + str(max_af)+ "' " + vcf_in 
            cmd_line += " | bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%TGT\t%GT\t%INFO/AF\n]' -i '" + 'GT = "alt"' + "'"


            print("Running this command")
            print(cmd_line)
            
            po = subprocess.Popen(cmd_line,
                                  shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

            stdout, stderr = po.communicate()

            po.wait()

            return_code = po.returncode

            if return_code != 0:
                raise Exception("Command line {} got return code {}.\nSTDOUT: {}\nSTDERR: {}".format(cmd_line, return_code, stdout, stderr))
            variants =  variants + [x.split('\t') for x in stdout.decode().split('\n')[:-1]]

    variants = [int(x[1]) for x in variants]

    print("Number of inherited controls: " + str(len(variants)))
    return variants
