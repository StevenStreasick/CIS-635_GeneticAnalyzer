# GeneticAnalyzer
## The Project 
For my project, I was particularly interested in analyzing genetics of some sorts. Some of my previous research involved analyzing genetics of high altitude populations such as the Tibetan and Andean population to other low land populations. My initial idea was to create an agent that would be able to identify people of high altitude populations. However, due to a combination of time and access limitations, this project quickly grew out of scope for this class. So an alternative project idea of creating a clustering algorithm based on the different populations was proposed and adopted. 

## Project Hypotheses
The core concept of this project relies on two hypotheses. The success of this project in identifying individuals origins would in turn prove these hypotheses, whereas a failure within this project would suggest one, or both, of the hypotheses being incorrect. 

Each population has undergone selective pressure upon their genetic data, causing the individuals to be ‘fitted’ to their environment. 
The number of variants and the kmer counting of the variants is enough to identify the population origins.

## The Data
The data was downloaded from the 1000 genomes project, under the hrc30 sampling. To save computer space, the variant data for all 2000+ individuals was downloaded. This was done using a Python call to bash utilizing the cURL language to iteratively download all 23 chromosome variant data. In an attempt to save memory, each file was downloaded, read from, and then deleted. 

## Data format
The downloaded files were under the .vcf.gz format, where VCF is considered to be a technique for storing massive medical data. The vcf.gz format is the compressed version of this data. This format is considered to be special, as it can be easily read and indexed using tabix to get full rows. This presents a new challenge given the structure of the data, where the column indicates the individual and the row indicates the variant position. Therefore gathering a full row would simply get all individuals variants for that specific position, indicating some sort of memory is required to store the previous row/rows of variant data to allow for kmer counting. To do this, I created a sliding window across the rows, creating new kmers as the sliding window shifts. Each time a new kmers is created, I update that individuals stored sparce matrix with the results from the kmer hashed using the hash vector. At the end, I'm left with m by n array of sparce matrix, where m is the number of unique tokens (kmers) found, and n is the number of individuals. 

## Hash Vectorizer
Researching several of the different count vectorizers, I chose to utilize the hash vectorizer because of it's utilization of the sparce matrix in addition to its fast speeds. The default hash function was utilized, as I had not a chance to research prior hash functions utilized in this field/application.

## Running the Application
I estimated this project running for approximately a day, as the variant data often downloaded in 20-30 minutes in my testing of the HPC. Given this, I doubled this to estimate the time to reach the file and utilize the counting. In practice, the program ran for 2 days before being stopped by the HPC because of the allocated time being exceeded. 

## Threats to Validity
One potential threat to validity is my utilization of kmer counting. Because of my implementation of it, the program treats the two combinations of variant pairs differently, and as such, could skew the model and results. Another potential threat to validity is the utilization of the hashing. Nothing ensures that the pairs are hashed the same way every time given a different occurrence pattern. This in turn could lead to sparse matrices being added that order their bins differently, leading to incorrect results 
