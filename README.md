# prosit_fasta2csv
Convert your fasta files to prosit compatible csv files.

### How to use this repo:
Clone the repo and copy&paste your fasta or fas file into the folder.
Change the name to "test_file.fasta", overwriting the provided test file. 
Start a terminal in the folder and run "python fasta2csv.py".

This will create a csv file in the prosit supported structure: (www.proteomicsdb.org/prosit)
where the charge states are set to 2 and collision engery of 28 is used.

To modify ce and precursor charge states, the fasta2csv function has to be run with
additional arguments. Use "python fasta2csv.py --help" to figure out how to do so.
example: python fasta2csv.py --ce 30 --cs 23
to set cillision energy to 30keV and randomly select charge states of 2 and 3. 
