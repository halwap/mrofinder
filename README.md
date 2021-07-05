# [mrofinder](https://github.com/halwap/mrofinder)

Automated pipeline for collecting various evidence for MRO-candidate proteins:
- homology based searches with:
  - hmm profiles of known MRO proteins, bidirectional blast with known mitochondrial proteins from [Mitominer](https://mitominer.mrc-mbu.cam.ac.uk/release-4.0/begin.do)
- structular analysis:
  - signal peptide prediction, identification of potential tail-anchored proteins and MBOMPs
- helpful in annotation searches:
  - blast against NCBI-nr, [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/)
      

## Requirements

- `Python >= 3.7`
- `biopython`
- `BLAST+`
- `hmmer`
- `targetp`
- `mitofates`
- `tmhmm`
- `psipred`

## Optional
- `InterProScan`


## Installation

```bash
git clone https://github.com/halwap/mrofinder
```
It is required to update the paths to all the external software and databases in mrofinder_paths.txt. It is possible, that one day proper conda package will be prepared :slightly_smiling_face:
Reminder: THE SOFTWARE IS PROVIDED "AS IS"


## Usage


#### Basic usage:
```bash
python /path/to/mrofinder_main.py -i input_proteins.fasta -o output_name.tsv (-w working_directory_name)
(--hmmer_results hmmer_results_dict/ | --hmmer_queries hmmer_queries_dict/)
(--q2db_blast protein2mitominer_blast.tab --db2q_blast mitominer2proein_blast.tab | --mitominer_fasta /path/to/mitominer.fasta --mitominer_db /path/to/mitominer_blast_db)
(| tee log.txt)
```


The results are presented in table with unnamed columns (for easier further handling).


mrofinder only works with protein sequences. Mitominer files are included in repository, and probably their paths should be defaulted to the repository #TODO

All intermediate files can be found in working directory (default working_directory_YYYY_MM_DD_HH_MM). 

The output columns are as following (if particular search was not conducted, corresponding column is absent):
protein_id, hmmer_hit, mitominer_hit, go_category_id, go_category_name, nr_hit_id, nr_hit_descrtiption, nr_hit_evalue,
N_of_software_detecting_peptide_signal, targetp_localisation, targetp_probability, mitofates_localisation, mitofates_probability, nommpred_localisation,
TA_or_not_TA, TA_helise_length, TA_helise_hydrophobicity, TA_before_charge, TA_end_charge, MBOMP_or_not_MBOMP



#### Advanced:

```bash
mrofinder --help
```


## License

mrofinder is released under an open-source MIT license 











