#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <fstream>
using namespace std;


float dist(float* v1, float* v2)
{
	return ((v1[0]-v2[0])*(v1[0]-v2[0])		+	(v1[1]-v2[1])*(v1[1]-v2[1])	+	 (v1[2]-v2[2])*(v1[2]-v2[2]));

}

int main(int argc, char const *argv[])
{
	
	cout<<"First parameter: orginal file. Second: Modified filename.\n ";	
	ifstream file(argv[1]);			
	
	ofstream file2(argv[2]);
	
	cout<<argv[1]<<" "<<argv[2]<<"\n";
	//getchar();
	char str[10];
	int dim;
	int noPoints;
	//file2<<"I am indside";
	file>>str;	dim = atoi(str);		
	file2<<str<<" ";
	file>>str; noPoints=atoi(str);
	cout<<"\n dim:"<<dim<<" noPoints:"<<noPoints;
	file2<<str<<"\n";
	std::vector<float*> inPoints = std::vector<float*>();

	for(int i=0;i<noPoints;i++)
	{	bool flag = false;
		float *xyz=(float*)malloc(dim*sizeof(float));

		for(int i=0;i<dim;i++)
		{
			
				file>>str;
				xyz[i] = atof(str);
				flag = true;
				file2<<str;
				if(i!=dim-1)
					file2<<" ";
			
			//cout<<xyz[i]<<" ";
		}
		file2<<"\n";
		if(flag == true)
		{
			inPoints.push_back(xyz);
			
		}
		
	}
		// all points sampled


	while(!file.eof())
	{
		file>>str;	//dimension
		file>>str;
		file2<<str<<" ";	//first coord
		int coord1 = atoi(str);
		file>>str;
		file2<<str<<" ";	//second coord
		int coord2 = atoi(str);
		file2<<dist(inPoints[coord1],inPoints[coord2])<<"\n";
		//file2<<"1\n";	//assuming weight is one


	}
	file.close();
	file2.close();
	return 0;
}