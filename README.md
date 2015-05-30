# facial-expression-recognition
This is a class project for the course of Data structure & algorithm. In this project, we were asked to implement a system to recognize human facial action units (AU) from facial images. Each image is represented by a predefined feature set (a vector containing 5,632 elements) and denoted by xi. Given a template image dataset Y = {y1, …, yN} (containing N = 138 images in this project). The recognition of AUs is realized by using a k-nearest neighbor-based method, by comparing the similarity between xi and each of the template images yj and finding the 10-nearest neighbors of xi. 
 
By using this similarity function, the most similar template image (the closest neighbor of xi) has the maximum similarity value.

The methodology I use to solve this problem includes three parts: 
1. Method to read input query and template files into one dimensional or two dimensional vectors for convenience of further computation. Six different routines to read file have been compared in terms of time efficiency: (a) using C API to read directly into a string and then read stringstream to double number vectors, (b) using C++ streams to read into string and convert to double number vectors, (c) using istreambuf_iterator to fast iteration out of stream buffers (files) and then convert to double numbers, (d) also using istream, but assign the size of the string first and then allocate the data into it, (e) fast copy files to another stream via operator << on their internal buffers and then handle the string to vector, (f) using getline function in istream and convert each line to a one dimensional vector and combine different lines to two dimensional vector.

2. Method to compute similarity value as well as optimize the process of computation. By brute force, it is just to compute the similarity value table of xi for the 138 different elements yj in template files.  New algorithm has been proposed by setting an upper bound to filter out unpromising elements and only compute similarity value for elements with higher upper bound. Computation efficiency for this algorithm is hard to analyze theoretically but improved computing time around 33% from experimental results.

3. Method to select the top 10 largest similarity value for each query file and sort them. A modified quicksort algorithm is applied to sort the 10 largest values from the 138 results: a pivot corresponding to the 10th largest value is being calculated using partial quicksort in the time O(n), and then another quicksort is taken to take care the first ten unsorted results, so the final efficiency comes out as O(n + klogk), k is the number of largest value, in this case 10.

Summary:
1.	Six different I/O algorithms were compared for reading large files into vectors
2.	An upper bound algorithm was designed to improve the time to compute cosine similarity function and improved 33% compared to brute force
3.	A partial quick sort algorithm is applied to sort the k largest number in a time O(n + klogk)

