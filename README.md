# test_submission
Submission for Biographica bioinformatics test

Q1. Open the task1 folder in the tasks folder. The data assembly_1.prot.fa and assembly_2.prot.fa are proteins from the same organism but different assemblies. Disambiguate the proteins based on sequence similarity (e.g. using relatively relaxed mapping of 95% sequence identity and 0.9 coverage factor). Therefore create a new fasta file containing all proteins from the original two fasta files but where each individual sequence only appears once and both IDs for that sequence included in the fasta record.


Q2. As can be seen in the file, many of the protein sequences for assembly 2 are standardized gene identifiers. Write a script to map the identifers, where possible, to uniprot. Therefore identify a protein in assembly_1.prot.fa which has hydrolase activity, acting on ester bonds.

## Notes:
- Please provide solutions in Python.
- We anticipate approximately 3hrs coding for the test above.
- Although not formally part of the coding test, we are also interested in enforcing data governance and metadata logging, as well as implementing pipeline versioning. Some of these considerations will go beyond the scope of the submission scripts. However, where possible, candidates should include in their submission code, comments detailing how they would scale their solution to millions of proteins. Where not possible we just expect candidates to have considered these and be willing to chat about them at the interview.
- For task 2 we don't at all intend for the task to take a prohibitive amount of time. If the building of a mapping script is prohibitively long, we can provide a script to automate ID mapping to uniprot identifiers. 

