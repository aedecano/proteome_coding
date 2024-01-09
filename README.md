# test_submission
Submission for Biographica bioinformatics test

**Q1.** Open the task1 folder in the tasks folder. The data assembly_1.prot.fa and assembly_2.prot.fa are proteins from the same organism but different assemblies. Disambiguate the proteins based on sequence similarity (e.g. using relatively relaxed mapping of 95% sequence identity and 0.9 coverage factor). Therefore create a new fasta file containing all proteins from the original two fasta files but where each individual sequence only appears once and both IDs for that sequence included in the fasta record.

**A1.** Inferring the pan proteome from protein assemblies is the optimal way to identify the shared and unique protein products between two or more organisms of the same species. To perform this task, run the custom code *get_pan_proteome.py* in this repo on the input protein assemblies. It mainly uses SeqIO and AlignIO from the biopython module for pairwise alignments.

Sample usage: Use the two files in test_assemblies as input for testing. 

```
python3 get_pan_proteome.py input_assembly1 input_assembly_2 output_consensus_pan_proteome.fasta
```
A sample output file from a run using bacterial proteomes (ecoli_consensus_panprot.fasta.gz) can be downloaded from this repo.

To scale the run to more than a pair to millions, use *multientry_panproteome.py* as follows. 

```
python3 multientry_panproteome.py /path/to/input_directory /path/to/output_pan_proteome.fasta
```

This script takes the input files in batches to reduce memory usage. 

**Q2.** As can be seen in the file, many of the protein sequences for assembly 2 are standardised gene identifiers. Write a script to map the identifers, where possible, to uniprot. Therefore identify a protein in assembly_1.prot.fa which has hydrolase activity, acting on ester bonds.

**A2.** This task involves annotating a protein assembly with ambiguous notation of proteins with the previously identified protein products documented in UniProt. To do this, run the code *annotate_uniprot.py* as below. 

```
python3 annotate_uniprot.py input_assembly_protein.fasta output.tsv
```

As for the target protein that functions as a hydrolase acting on ester bonds: look for the myrosinase or myrosinase-like protein in the output.tsv.

NB: I used the assemblies provided by Dominic as input but the runs took a very long time and were eventually killed in our cluster.


## Notes:
- Please provide solutions in Python.
- We anticipate approximately 3hrs coding for the test above.
- Although not formally part of the coding test, we are also interested in enforcing data governance and metadata logging, as well as implementing pipeline versioning. Some of these considerations will go beyond the scope of the submission scripts. However, where possible, candidates should include in their submission code, comments detailing how they would scale their solution to millions of proteins. Where not possible we just expect candidates to have considered these and be willing to chat about them at the interview.
- For task 2 we don't at all intend for the task to take a prohibitive amount of time. If the building of a mapping script is prohibitively long, we can provide a script to automate ID mapping to uniprot identifiers.
