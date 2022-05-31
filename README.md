UVC-delins is a very accurate and reasonably fast delins (including multiple-nucleotide variant (MNV) and complex insertion-deletion (InDel)) variant caller based on the VCF file generated by UVC (which is available at https://github.com/genetronhealth/uvc). 
Please note that UVC already stored each haplotype and its information during the process of calling simple small variants. 
In other words, the step "store each haplotype and its information during the process of calling simple small variants" was already done in UVC. 
Therefore, the VCF file generated by UVC already contains haplotype information which is used by UVC-delins to call deletion-insertion (delins) variants. 

In fact, UVC-delins can process any VCF file with variant records containing parenthesis-grouped haplotype IDs (for example, 'Hap=(haplotype1)(haplotypep2)' without single quotes) in the VCF tags/fields of these records, meaning that this variant record is covered by haplotype1 and haplotype2. 
For more information about the forms of haplotype IDs used by UVC, please see the description for the bHap VCF tag in UVC. 

Although UVC-delins can be used with other SNV/InDel variant callers, we highly recommend to use it with UVC due to its seamless integration with UVC. 

# How to install

The installation requirements for UVC-delins are the same as these for UVC (for the requirements for UVC, please check https://github.com/genetronhealth/uvc)

In sum, UVC-delins requires BASH 4.0+ (4.0 is the minimum version required) and a compiler that supports the C++14 standard. The Makefile in this directory compiles with g++, but the Makefile can be easily modified to use another compiler instead of g++ (for example, clang). To install from scratch, please run: (./install-dependencies.sh && make clean && make all -j4 && make deploy).

In total, the installation for UVC-delins should take about 5 minutes.

# How to use

The script uvcvcf-raw2delins-all.sh in the bin directory is the main executable that generates all VCF files related to the calling of delins variants.
Run bin/uvcvcf-raw2delins-all.sh without any command-line argument will display its usage help.
The usage help for uvcvcf-raw2delins-all.sh refers to the executable uvcvcf-raw2delins, which performs the actual variant calling.
The binary executable uvcvcf-raw2delins performs the actual calling of delins variants by combining SNV(s) and InDel(s) that are near each other from the same haplotype, where the SNV(s) and InDel(s) that are not from the same haplotype are still kept.  
The script uvcvcf-raw2delins-all.sh simply wraps around the binary executable uvcvcf-raw2delins.

In sum, given an input BAM file, the command to generate VCF with delins variants (which include MNVs and complex InDels) along with SNVs and simple InDels are as follows: 

UVC_INSTALL_DIRECTORY/uvc1 -f INDEXED_REFERENCE_FASTA -o SNV_INDEL_VCF_GZ INPUT_BAM && UVC_DELINS_INSTALL_DIRECTORY/uvcvcf-raw2delins-all.sh INDEXED_REFERENCE_FASTA SNV_INDEL_VCF_GZ DELINS_VARIANT_VCF_PREFIX

# What to report if a runtime error arises

All bug reports, feature requests, and ideas for improvement are welcome (although not all of them may be addressed in time)!

# Other things

For more information, please check the wiki.

# References

## Publication

TODO

## Patent

TODO

