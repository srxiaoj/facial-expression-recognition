/*This program is used to select the most 10 objects in the template that
resemble the query file   ver 1.3 */ 
//header files #include <iostream>
#include <iostream>
#include <fstream> 
#include <string> 
#include <vector> 
#include <time.h>
#include <math.h> 
#include <algorithm>
#include <sstream>
#include <cerrno>
using std::vector;      
using std::cout; 
using std::string;      
using std::ofstream; 
using std::ifstream;    
using std::ios; 
using std::endl;
using std::swap;

//define constant numbers
int const numOfElement = 5632;
int const numOfVectors = 138;
int const numberOfSort = 10;

//prototypes
struct Svalue;
vector <Svalue> getSvalue(vector<double> &aVector);
string fileToString(string filename);
vector <double> generateOneDVector(string filename);
vector < vector<double> > generateTwoDVector(string filename);
void insertionSort(vector<Svalue> & aSvalue, int start, int end);
void selectionSort(vector<Svalue> & aSvalue, int start, int end, int k);
//Svalue quickSort(vector <Svalue> &aSvalue, int start, int end, int k);
void quickSort(vector <Svalue> &aSvalue, int l, int r);
int partialQuickSort(vector <Svalue> &aSvalue, int start, int end, int k);
int partition(vector <Svalue> &aSvalue, int start, int end);
void printVector(vector<double> aVector, int size);
void printStruct(vector <Svalue> aSvalue);
void outputOneDVector(vector <double> aOneDVector, string filename);
void outputTwoDVector(vector <vector <double> > aTwoDVector);
void outputSTable(vector <Svalue> aSvalue, ofstream &outputFile);
void outputSValue(vector <Svalue> aSvalue, string queryFileName, ofstream & outputFileName);
void outputSValue(vector <int> aSvalue, string queryFileName, ofstream &outputFile);
double getSumXSquare (vector <double> &xVector);
vector <double> getSumYSquare(vector <vector <double> > &yTwoDVector);
vector <double> getSumXY(vector <double> &xVector, vector <vector <double> > &yTwoDVector);
vector <double> getSimilarity(double &sumXSquare, vector <double> &sumXYVector, vector <double> &sumYSquareVector);
vector <string> getQueryFileNames(string auIndex);
vector <string> getTemplateFileNames();
vector <Svalue> getSVector(vector <double> &queryVector, vector<vector <double> > &templateVector, vector <double> &sumYSquare);
vector <Svalue> getTopValue(vector <Svalue> &svalue, int k);

int main()
{
	time_t startTime = time(NULL);//recording start time
	
	vector <string> AU01_FileName = getQueryFileNames("_AU01");
	vector <string> AU12_FileName = getQueryFileNames("_AU12");
	vector <string> AU17_FileName = getQueryFileNames("_AU17");
	vector <string> templateFileName = getTemplateFileNames();	
	ofstream outputFile("Similarity_Table.dat");
	outputFile.close();
	outputFile.open("Similarity_Table.dat", ios::out|ios::app);
	cout << "start" << endl;
	double inputTime = 0;
	double computeTime = 0;
	double sortTime = 0;
	for(int i = 0; i < 47; i++)
	{
		double beginInput = time(NULL);
		//read data from query file, store as a vector with 5632 row
		vector <double> queryVector = generateOneDVector(AU01_FileName[i]);
		//read data from template file, store as two vectors with 138 column and 5632 row
		vector<vector <double> > templateVector = generateTwoDVector(templateFileName[i]);
		inputTime += time(NULL) - beginInput;
		//compute sum(yi,k*yi,k)
		double beginCompute = time(NULL);
		vector <double> sumYSquare = getSumYSquare(templateVector);
		vector <Svalue> svalue = getSVector(queryVector, templateVector, sumYSquare);
		computeTime += time(NULL) - beginCompute;
		double beginSort = time(NULL);
		svalue = getTopValue(svalue, numberOfSort);
		outputSValue(svalue, AU01_FileName[i], outputFile);
		sortTime += time(NULL) - beginSort;
		svalue.clear();

		beginInput = time(NULL);
		queryVector = generateOneDVector(AU12_FileName[i]);
		inputTime += time(NULL) - beginInput;
		beginCompute = time(NULL);
		svalue = getSVector(queryVector, templateVector, sumYSquare);
		computeTime += time(NULL) - beginCompute;
		beginSort = time(NULL);
		svalue = getTopValue(svalue, numberOfSort);
		outputSValue(svalue, AU12_FileName[i], outputFile);
		sortTime += time(NULL) - beginSort;
		svalue.clear();

		beginInput = time(NULL);
		queryVector = generateOneDVector(AU17_FileName[i]);
		inputTime += time(NULL) - beginInput;
		beginCompute = time(NULL);
		svalue = getSVector(queryVector, templateVector, sumYSquare);
		computeTime += time(NULL) - beginCompute;
		beginSort = time(NULL);
		svalue = getTopValue(svalue, numberOfSort);
		outputSValue(svalue, AU17_FileName[i], outputFile);
		sortTime += time(NULL) - beginSort;
		svalue.clear();
	}

	//double readTime = time(NULL) - startTime;
	//cout << "Read time is: " << readTime << " seconds" << endl;
	
	//output time taken
	double completeTime = time(NULL) - startTime;
	int second = fmod(completeTime, 60);
	int min = (completeTime - second)/60;
	cout << "total time: " << min << " minutes " << second << " seconds" << endl;
	cout << "Input time: " << inputTime << ", Compute time: " << computeTime << ", Sort time: " << sortTime << endl;
	return 0;
}

/*read data from query file, store as a vector with 5632 row*/
vector <double> generateOneDVector(string filename)
{
	/*string aString = fileToString(filename);
	std::stringstream stream(aString);
	vector <double> aOneDVector;
	double aValue;
	while (stream >> aValue)
		{
			//cout << "The value is: " << aValue << endl;
			aOneDVector.push_back(aValue);
			if (stream.peek() == ' ')
				stream.ignore();
		}
		return aOneDVector;*/

	ifstream aFile (filename.c_str());
	vector <double> aOneDVector;
	double aValue;
	if (aFile.is_open())
	{
		while (aFile >> aValue)
		{
			//cout << "The value is: " << aValue << endl;
			aOneDVector.push_back(aValue);
			if (aFile.peek() == ' ')
				aFile.ignore();
		}
		aFile.close();
	}
	else 
		cout << "Unable to open file" << endl;
	return aOneDVector;

}
/*convert file to string*/
string fileToString(string filename)
{
	/*std::FILE *fp = std::fopen(filename.c_str(), "rb");
  	if (fp)
  	{
    	std::string contents;
    	std::fseek(fp, 0, SEEK_END);
    	contents.resize(std::ftell(fp));
    	std::rewind(fp);
    	std::fread(&contents[0], 1, contents.size(), fp);
    	std::fclose(fp);
    return(contents);
    }
    throw(errno);*/
   
    std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
  	if (in)
  	{
    	std::string contents;
    	in.seekg(0, std::ios::end);
    	contents.resize(in.tellg());
    	in.seekg(0, std::ios::beg);
    	in.read(&contents[0], contents.size());
    	in.close();
    	return(contents);
    }
    throw(errno);

    /*std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
    if (in)
   	{
    	std::ostringstream contents;
    	contents << in.rdbuf();
    	in.close();
    	return(contents.str());
    }
    throw(errno);*/
}
/*read data from template file, store as two vectors with 138 column and 5632 row*/
vector < vector<double> > generateTwoDVector(string filename)
{
	/*string aString = fileToString(filename);
	std::stringstream stream(aString);
	vector <vector <double> > aTwoDVector;
	double aValue;
	for (int row = 0; row < numOfVectors; row++)
	{
		vector <double> aOneDVector;
		aOneDVector.clear();
		for(int column = 0; column < numOfElement; column++)
		{
			stream >> aValue;
			aOneDVector.push_back(aValue);
		}
		aTwoDVector.push_back(aOneDVector);
	}
	return aTwoDVector;*/
	
	string aString = fileToString(filename);
	std::stringstream stream(aString);
	vector <vector <double> > aTwoDVector;
	double aValue;
	while(stream >> aValue)
	{
				vector <double> aOneDVector;
				aOneDVector.clear();
				while(stream.peek() != '\n')
				{
					aOneDVector.push_back(aValue);
					//cout << "the value is: " << aValue << endl;
					if (stream.peek() == ' ')
					{
						stream.ignore();
						stream >> aValue;
					}
				}
				aOneDVector.push_back(aValue);
				aTwoDVector.push_back(aOneDVector);
	}
	return aTwoDVector;

	/*ifstream aFile (filename.c_str());
	vector <vector <double> > aTwoDVector;
	double aValue;
	if(aFile.is_open())
	{
			while(aFile >> aValue)
			{
				vector <double> aOneDVector;
				aOneDVector.clear();
				while(aFile.peek() != '\n')
				{
					aOneDVector.push_back(aValue);
					//cout << "the value is: " << aValue << endl;
					if (aFile.peek() == ' ')
					{
						aFile.ignore();
						aFile >> aValue;
					}
				}
				aOneDVector.push_back(aValue);
				aTwoDVector.push_back(aOneDVector);	
			}
		aFile.close();
	}
	else
		cout << "Unable to open file" << endl; 
	return aTwoDVector;*/

	/*
	string aString;
	vector<vector <double> > aTwoDVector;
	size_t posOfSpace = 0;
	if(aFile.is_open())
	{
		int row = 0;
		while(getline(aFile, aString))
		{
			vector<double> oneDVector;//save each row in a oneD vector
			int column = 0;
			while(!aString.empty())
			{
				oneDVector.push_back(stod(aString));
				posOfSpace = aString.find(" ");
				if(posOfSpace<numOfElement)
				{
					aString = aString.substr(posOfSpace+1);
					column++;
				}
				else aString = "";
			}
			aTwoDVector.push_back(oneDVector);
		}
		aFile.close();
	}
	else cout << "Unable to open file" << endl;
	return aTwoDVector;
	*/
}
/* a structure type contains value and index */
struct Svalue
{
	double value;
	int index;
};
/*convert a vector to a vector contains value and index */
vector <Svalue> getSvalue(vector<double> &aVector)
{
	vector <Svalue> aSvalue;
	for(int i = 0; i<aVector.size(); i++)
	{
		Svalue *newSvalue = new Svalue;
		newSvalue -> value = aVector[i];
		newSvalue -> index = i;
		aSvalue.push_back(*newSvalue);
	}
	return aSvalue;
}
/*Partial QuickSort: find the kth largest element in the vector*/
int partialQuickSort(vector <Svalue> &aSvalue, int start, int end, int k)
{
	int pivot;
	if(end > start)
	{
		pivot = partition(aSvalue, start, end);
		if(pivot == k-1)
		{
			//cout << "pivot is: " << pivot << " ";
			return aSvalue[pivot].index;
		}
		else if(pivot > k - 1)
		{
			partialQuickSort(aSvalue, start, pivot, k);
		}
		else
		{
			partialQuickSort(aSvalue, pivot + 1, end, k);
		}
	}
	else return aSvalue[end].index;
}
/* quicksort*/
void quickSort(vector <Svalue> &aSvalue, int l, int r)
{
	int pivot;
	if (r > l)
	{
		pivot = partition(aSvalue, l, r);
		quickSort(aSvalue, l, pivot - 1);
		quickSort(aSvalue, pivot + 1, r);
	}
}
/*partition*/
int partition(vector <Svalue> &aSvalue, int left, int right)
{
	int p = left;
	double pValue = aSvalue[left].value;
	for (int i = left+1; i <= right; i++)
	{
		if(aSvalue[i].value > pValue)
		{
			p++;
			swap(aSvalue[i].value, aSvalue[p].value);
			swap(aSvalue[i].index, aSvalue[p].index);
		}
	}
	swap(aSvalue[p].value, aSvalue[left].value);
	swap(aSvalue[p].index, aSvalue[left].index);
	return p;
	/*int i = left, j = right;
	cout << "i:" << i << " j: "<< j << endl;
	while (j > i)
	{
		while (aSvalue[i].value < aSvalue[p].value)
		{
			i++;
		}
		while (aSvalue[j].value > aSvalue[p].value)
		{
			j--;
		}
		cout << "swap" << " ";
		swap(aSvalue[i].value, aSvalue[j].value);
		swap(aSvalue[i].index, aSvalue[j].index);
	}
	swap(aSvalue[i].value, aSvalue[j].value);
	swap(aSvalue[i].index, aSvalue[j].index);
	p = j;
	return p;*/
}

/*insertion sort*/
void insertionSort(vector<Svalue> & aSvalue, int start, int end)
{
    int i, j;
    for(i = start+1; i < end; i++)    // Start with 1 (not 0)
    {
        j = i;
        while (j > 0 && aSvalue[j].value > aSvalue[j-1].value)
        {
        	swap(aSvalue[j].value, aSvalue[j-1].value);
        	swap(aSvalue[j].index, aSvalue[j-1].index);
        	j--;
        }
    }
    return;
}
/*Selection sort:  Sort a vector, using the selection sort algorithm, on the order of O(10*138) */
void selectionSort(vector<Svalue> & aSvalue, int start, int end, int k)
{
    for (int i = start; i < k; i++)//only sort the first ten value
    {
        int min = i;
        for (int j = i+1; j < end; j++)
        {
			if (aSvalue[j].value < aSvalue[min].value)
            {
                min = j;
            }
        }

        if (min != i)
        {
            swap(aSvalue[i].value, aSvalue[min].value);
			swap(aSvalue[i].index, aSvalue[min].index);
        }
		//cout << "min index: " << aSvalue[min].index << endl;
    }
}

/*input: query file name, tempalte file name, output: sorted similarity value */
/*vector <int> getSVector(string queryFileName, string templateFileName)
{
	//read data from query file, store as a vector with 5632 row
	vector <double> queryVector = generateOneDVector(queryFileName);
	//outputOneDVector(vector_001_AU01);
	//read data from template file, store as two vectors with 138 column and 5632 row
	vector<vector <double> > templateVector = generateTwoDVector(templateFileName);
	
	//test if template date is stored correctly
	//outputTwoDVector(twoDVector_001_template);
	
	//compute sum(xk*xk)
	double sumXSquare = getSumXSquare(queryVector);
	//compute sum(yi,k*yi,k)
	vector <double> sumYSquare = getSumYSquare(templateVector);
	//compute sum(xk*yi,k)
	vector <double> sumXY = getSumXY(queryVector, templateVector);
	//compute similarity values
	vector <double> s = getSimilarity(sumXSquare, sumXY, sumYSquare);
	//cout << "Complete similarity function" << endl;
	//test if s is calculated correctly
	//OutputOneDVector(s);

	//presort similarity values
	vector <Svalue> svalue = getSvalue(s);
	//quickSort(svalue,0,svalue.size(), 2);
	vector <int> sortedS;
	sortedS.clear();
	int value = quickSort(svalue, 0, svalue.size(), numberOfSort);
	selection(svalue, 0, numberOfSort);
	return sortedS;
}
*/
/*input: query file name, tempalte file name, output: sorted similarity value */
vector <Svalue> getSVector(vector <double> &queryVector, vector<vector <double> > &templateVector, vector <double> &sumYSquare)
{	
	//compute sum(xk*xk)
	double sumXSquare = getSumXSquare(queryVector);
	//compute sum(xk*yi,k)
	vector <double> sumXY = getSumXY(queryVector, templateVector);
	//outputOneDVector(sumYSquare, "sumYSquare.dat");
	//compute similarity values
	vector <double> s = getSimilarity(sumXSquare, sumXY, sumYSquare);
	//cout << "Complete similarity function" << endl;
	//test if s is calculated correctly
	//OutputOneDVector(s);

	//presort similarity values
	vector <Svalue> svalue = getSvalue(s);
	return svalue;
}
/*get the top ten result*/
vector <Svalue> getTopValue(vector <Svalue> &svalue, int k)
{
	//quickSort(svalue,0,svalue.size(), 2);
	int value = partialQuickSort(svalue, 0, svalue.size(), k);
	quickSort(svalue, 0, k);
	//selectionSort(svalue, 0, svalue.size());
	return svalue;
}
/*compute sum(xk*xk)*/
double getSumXSquare (vector <double> &xVector)
{
	double sumXSquare = 0.0;
	for(int k = 0; k < numOfElement; k++)
		{
			sumXSquare += xVector[k]*xVector[k];
		}
	sumXSquare = sqrt(sumXSquare);
	//cout << "sum X: " << sumXSquare << endl;
	return sumXSquare;
}
/*compute sum(yi,k*yi,k) */
vector <double> getSumYSquare(vector <vector <double> > &yTwoDVector)
{
	vector <double> sumYSquare;
	for(int i = 0; i<numOfVectors; i++)
	{
		double count = 0.0;
		for(int k = 0; k < numOfElement; k++)
		{
			count += yTwoDVector[i][k]*yTwoDVector[i][k];
		}
		count = sqrt(count);
		sumYSquare.push_back(count);
	}
	return sumYSquare;
}
/*compute sum(xk*yi,k) */
vector <double> getSumXY(vector <double> &xVector, vector <vector <double> > &yTwoDVector)
{
	vector <double> sumXY;
	for(int i = 0; i<numOfVectors; i++)
	{
		double count = 0.0;
		for(int k = 0; k < numOfElement; k++)
		{
			count += xVector[k]*yTwoDVector[i][k];
		}
		sumXY.push_back(count);
	}
	//outputOneDVector(sumXY, "sumXY.dat");
	return sumXY;
}
/*compute similarity function*/
vector <double> getSimilarity(double &sumXSquare, vector <double> &sumXYVector, vector <double> &sumYSquareVector)
{
	vector <double> s;
	for(int i = 0; i<numOfVectors; i++)
	{
		s.push_back(sumXYVector[i]/(sumXSquare*sumYSquareVector[i]));
	}
	return s;
}
/*output oneD vector*/
void outputOneDVector(vector <double> aOneDVector, string filename)
{	
	ofstream oneDFile(filename.c_str());
	for (int i = 0; i < aOneDVector.size(); i++)
	{
		oneDFile << aOneDVector[i] << " ";
	}
	oneDFile.close();
}
/*output twoD vector*/
void outputTwoDVector(vector <vector <double> > aTwoDVector)
{
	ofstream testTemplate("TwoDVector.dat");
	for(int i = 0; i<aTwoDVector.size(); i++)
	{
		for(int j = 0; j<aTwoDVector[i].size(); j++)
		{
			testTemplate << aTwoDVector[i][j] << " ";
		}
		testTemplate << endl;
	}
	testTemplate.close();
}
/*output s value with index*/
void outputSTable(vector <Svalue> aSvalue, ofstream &outputFile)
{
	for (int i = 0; i < aSvalue.size(); i++)
	{
		outputFile << aSvalue[i].value << ",";
		outputFile << aSvalue[i].index << " ";
	}
	outputFile << '\n';
}

/*output s value table*/
void outputSValue(vector <Svalue> aSvalue, string queryFileName, ofstream &outputFile)
{
	outputFile << queryFileName << ": ";
	for (int i = 0; i < numberOfSort; i++)
	{
		outputFile << aSvalue[i].value << " ";
		outputFile << aSvalue[i].index << ", ";
	}
	outputFile << '\n';
}
/*output s value table*/
void outputSValue(vector <int> aSvalue, string queryFileName, ofstream &outputFile)
{
	outputFile << queryFileName << ": ";
	for (int i = 0; i < numberOfSort; i++)
	{
		outputFile << aSvalue[i] << " ";
	}
	outputFile << '\n';
}
/* PrintVector:  Print contents of an array.*/
void printVector(vector<double> aVector, int size)
{
	for (int i = 0; i < size; i++)
	{
		cout << aVector[i] << " ";
	}
	cout << '\n';
}
/* Print contents of a s value structure.*/
void printStruct(vector <Svalue> aSvalue)
{
	//print out the s value with index
	for(int i = 0; i<aSvalue.size(); i++)
	{
		cout << "(" << aSvalue[i].value << ", ";
		cout << aSvalue[i].index << ") ";
	}
}

/*get a vector to save all the query file name*/
vector <string> getQueryFileNames(string auIndex)
{
	vector <string> fileName;
	fileName.clear();
	string noIndex = ("_query.dat");
	for (int i = 1; i < 48; i++)
	{
		char indexString[10] = {""};
		sprintf(indexString, "%03d", i);
		string newFileName = string(indexString) + auIndex + noIndex;
		fileName.push_back(newFileName);
		//cout << "File name is: " << fileName[i-1] << endl;
	}
	return fileName;
}
/*get a vector to save all the template file name*/
vector <string> getTemplateFileNames()
{
	vector <string> fileName;
	fileName.clear();
	string noIndex = ("_template.dat");
	for (int i = 1; i < 48; i++)
	{
		char indexString[10] = {""};
		sprintf(indexString, "%03d", i);
		string newFileName = string(indexString) + noIndex;
		fileName.push_back(newFileName);
		//cout << "File name is: " << fileName[i-1] << endl;
	}
	return fileName;
}
