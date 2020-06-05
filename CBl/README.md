#Conway and Bromage succint data structures for assembling large genomnes

##Current content

### encode function
**encode** allow us to translate a nucleotide sequence into an int thanks to a given code :
A = 0
C = 1
G = 2
T = 3
We use an hash to create and return unique int. This will help us to complete the sd_vector
The function returns an int that matches with the location of a one in the future sd_vector


###fromFileToSdVector function
The goal of this function is to get sort nucleotide sequences which are store in a tubular file to transcribe them
into a sd_vector (from the sdsl library).
We call **encode** in **fromFileToSdVector** and we create the sd_vector from the location of ones information
At the end, we return an sd_vector which contain 0 and 1. Ones correspond to the location of elements which are in the file (after encode to understand where the location is).
