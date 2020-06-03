# Fake kmer generators

* generate_counts.py: This script can be used to generate fake kmers counts in a tubular file.
This is the simpliest generator that can be created.
All the counts are 1 and reverse complement are not canonized.
The kmers are **not** sorted.

Examples:

```bash
  # Generate 50 kmers of size 13 and output them on stdout
  python3 generate_counts.py -k 13 -n 50

  # Generate sorted file of 100 31-mers
  python3 generate_counts.py -k 31 -n 100 | sort > sorted_kmers.txt
```