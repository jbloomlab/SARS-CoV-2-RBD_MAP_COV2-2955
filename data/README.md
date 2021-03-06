# Input data

This directory contains input data used for the analysis and configuration for some of the analyses:

 - [barcode_runs.csv](barcode_runs.csv): Illumina sequencing runs for barcode counting.
   The *sample* column should be the hyphen separated concatenation of *experiment*, *antibody*, *concentration*, and *selection*.
   The **combination** of the *library* and *sample* columns should be unique.
   The *frac_escape* column gives the overall fraction of the library that escaped antibody.
   The *R1* FASTQ files should be semicolon (`; `) separated list of glob patterns.

 - [wildtype_sequence.fasta](wildtype_sequence.fasta): wildtype sequence of mutagenized gene.

 - [./pdbs/](pdbs): files downloaded from the [Protein Data Bank](https://www.rcsb.org/) for structural analyses.

 - [escape_profiles_config.yaml](escape_profiles_config.yaml): Information on how to plot the escape profiles; manually edit this to alter their plotting.

 - [output_pdbs_config.yaml](output_pdbs_config.yaml): Information on how to output structural mappings of escape.

 - [site_color_schemes.csv](site_color_schemes.csv): Schemes for how to color sites (can be used in escape profiles). Here are details on these schemes.

   - The `subdomain` scheme colors sites orange if they directly contact ACE2 (within 4 angstroms in PDB 6m0j) in the SARS-CoV-2 structure (residues 417, 446, 449, 453, 455, 456, 475, 486, 487, 489, 493, 496, 498, 500, 501, 502, 505), blue if they are in the receptor binding motif (RBM, residue 437 to 508, inclusive), and green if they are in the core RBD domain (all other sites). These definitions match those used in [Starr et al (2020)](https://www.cell.com/cell/fulltext/S0092-8674(20)31003-5).

 -  [spikeprot1030.fasta.tar.tar.xz](spikeprot1030.fasta.tar.tar.xz): A compressed FASTA file with all spike sequences in GISAID. Downloaded on Nov-03-2021 as follows: logged into GISAID, went to the *EpiCoV* page, clicked on *Downloads*, clicked on *spikeprot1030*, which yielded the file `spikeprot1030.fasta.tar.tar.xz` which was then copied here. We acknowledge all contributors to the GISAID database.

 - [./RBD_sites.csv](RBD_sites.csv): Useful numeric and functional annotations of RBD sites for referencing in various analyses.

 - [./RBDs_aligned.fasta](RBDs_aligned.fasta): alignment of sarbecovirus RBDs
