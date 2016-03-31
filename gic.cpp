/*
(c) 2012 Fengtao Fan
*/
#include "ANNSearchSampling.h"
#include "SimpleGraph.h"
#include "PointSet.h"
#include "GIComplex.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <boost/program_options.hpp>
#define ROUNDOFF 100
#include <float.h>
float rounding(float x)	{return ceil(x*ROUNDOFF)/ROUNDOFF;}
float minuss(float x, float y){ return abs(rounding(rounding(x)-rounding(y)));}
////////////////


// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int minDistance(float dist[], int V, bool sptSet[])
{
   // Initialize min value
   float min = FLT_MAX;
   int min_index;
 
   for (int v = 0; v < V; v++)
     if (sptSet[v] == false && dist[v] <= min)
         {min = dist[v]; min_index = v;}
 
   return min_index;
}
 
// A utility function to print the constructed distance array
int printSolution(int dist[], int V, int n)
{
   printf("Vertex   Distance from Source\n");
   for (int i = 0; i < V; i++)
      printf("%d \t\t %d\n", i, dist[i]);
}
 
// Funtion that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
int dijkstra(float *graph[], int V, int src, float* dist, std::vector<int> SubptInd)
{
	
    //static int dist[V];     // The output array.  dist[i] will hold the shortest
                      // distance from src to i
 
     bool sptSet[V]; // sptSet[i] will true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized
 
     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V; i++)
        {
        	//if(i!=src)
        	dist[i] = FLT_MAX;
        	sptSet[i] = false;
        }
 	
     // Distance of source vertex from itself is always 0
     dist[src] = 0;
 
     // Find shortest path for all vertices
     for (int count = 0; count < V-1; count++)
     {
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in first iteration.
     	//cout<<"inside dijkstra1";
		//getchar();
       int u = minDistance(dist, V, sptSet);
 		//cout<<"inside dijkstra2";
		//getchar();
       // Mark the picked vertex as processed
        sptSet[u] = true;
        if(u!=src && SubptInd[u] != 0)	//closes subsampled point
        	{if(SubptInd[src]!=0 ) 
        		{
        			cout<<"Error";	exit(0);	
        		}
        		return SubptInd[u];
        	}
 		//cout<<"inside dijkstra3";
		//																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																										getchar();
       // Update dist value of the adjacent vertices of the picked vertex.
       for (int v = 0; v < V; v++)
 		{
         // Update dist[v] only if is not in sptSet, there is an edge from 
         // u to v, and total weight of path from src to  v through u is 
         // smaller than current value of dist[v]
       	//if(count == 2 && v>6700) {	cout<<"inside dijkstraloop. u:"<<u<<" "<<v<<"\n";		getchar();}
         if (!sptSet[v] && graph[u][v] && dist[u] != FLT_MAX 
                                       && dist[u]+graph[u][v] < dist[v])
            dist[v] = dist[u] + graph[u][v];
    	}
     }

 return -1;
     // print the constructed distance array
   //  printSolution(dist, V);
}


////////////////
 int getSubPointIndex(const char* filename1,const char* filename2, std::map<int, int> &SubptInd)
{
	std::vector<float*> inPoints = std::vector<float*>();
	char str[10];
	int dim;	
	int noPoints;
	//std::map<int,int> SubptInd;

	/////////////////////////// gets all points from datapoints2///////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	cout<<" \nfrom get sub points.";
	ifstream file(filename1);			//datapoints2
	file>>str;	dim = atoi(str);		
	file>>str; noPoints=atoi(str);
	
	for(int i=0;i<noPoints;i++)
	{	bool flag = false;
		float *xyz=(float*)malloc(dim*sizeof(float));

		for(int i=0;i<dim;i++)
		{
			
				file>>str;
				xyz[i] = atof(str);
				flag = true;
			
			//cout<<xyz[i]<<" ";
		}
		if(flag == true)
		{
			inPoints.push_back(xyz);
			
		}
		
	}

	//inPoints.erase(inPoints.begin()+inPoints.size()-1);
	
	//for(int i=0;i<inPoints.size();i++)
	//	cout<<(inPoints[i])[0]<<" "<<(inPoints[i])[1]<<" "<<(inPoints[i])[2]<<"\n";
	//exit(0);
	
	file.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////		get all subsampled PointSet//////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	std::vector<float*> subPoints = std::vector<float*>();
	char stri[10];
	int dim2;	
	////////////////////////	Sizes
	int freqP;	int freqTot = inPoints.size();
	ifstream file1(filename2);
	file1>>stri;	dim2 = atoi(stri);
	file1>>stri; freqP = atoi(stri);
	cout<<" freq of subsampled points: "<<freqP<<" ";
	for(int j=0;j<freqP;j++)
	{
		float *xyz=(float*)malloc(dim*sizeof(float));
		for(int i=0;i<dim2;i++)
		{
			
				file1>>str;
				xyz[i] = atof(str);
			
			
		}
		
			subPoints.push_back(xyz);
		
	}
		
	//for(int i=0;i<subPoints.size();i++)
	//		cout<<(subPoints[i])[0]<<" "<<(subPoints[i])[1]<<" "<<(subPoints[i])[2]<<"\n";
	
	file1.close();
	//exit(0);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////


	for(int i=0;i<freqTot;i++)			//Initialise
			SubptInd[i]=0;
	int j=0;
	int countre = 0;
	cout<<" inside submap ";
	 for(int i=0;i<freqP;i++)
	 {
	// 	//subPoints[i]
	 	j=0;
	 	while(true)
	 		{
	 			float x = (subPoints[i][0]-inPoints[j][0]);
	 			float y = (subPoints[i][1]-inPoints[j][1]);
	 			float z = (subPoints[i][2]-inPoints[j][2]);
	 			if(x<0)	x = -x;	if(y<0) y = -y; if(z<0)	z =-z;
	 			if( x <=0.00001 && y <=0.00001 &&  z <= 0.00001)
	 				{
	 					//if(i==0)
	 					if(i==0)
	 						SubptInd[j] = freqTot+10;
	 							//SubptInd.insert ( std::pair<int,int>(j,i) );
	 							else SubptInd[j] = i;
	 							//->first = 

	 							//SubptInd[j]=freqTot+10;
	 					
	 					//else
	 					//	SubptInd[j]=i;
	 					countre++;
	 					break;
	 				}
	 			j++;
	 			if(j==freqTot)
	 				{
	 					cout<<"No match found "<<rounding(subPoints[i][0])<<" "<<rounding(subPoints[i][1])<<" "<<rounding(subPoints[i][2])<<"\n";;
	 					exit(0);
	 				}
	 		}
	 		//cout<<"matched i:"<<i<<" j:"<<j<<" ";
	 		//cout<<(inPoints[j])[0]<<" "<<(inPoints[j])[1]<<" "<<(inPoints[j])[2]<<"		WITH 	"; 	
			//cout<<subPoints[i][0]<<" "<<subPoints[i][1]<<" "<<subPoints[i][2]<<"\n";
			
	 	


	 }

cout<<countre<<". out of sub map\n";
return freqTot;
//vector<float*>().swap(inPoints);
//vector<float*>().swap(subPoints);
//return SubptInd;
	//vector<float*>().swap(inPoints);
	//vector<float*>().swap(subPoints);
}

/////////////////////////////////////////////////

//int findmin(int first[], int noPoints,std::map<int, int> SubptInd )
//{/
//	int min = INT_MAX;
//	int pos = -1;
//	for(int i=0;i<noPoints;i++)
//		if(first[i]<min && SubptInd[i]!=0)
//			{
//				min = first[i];
//				pos = i;
//			}
//return pos;
//}

/////////////////////////////////////////////////
void inputGraph(const char* filename1,const char* filename2, std::map<int, int> &SubPointIndex, std::vector<int> &ColorMapping )
{
	ifstream file(filename1);
	int dim;	
	int noPoints;
	char str[10];
	
	std::vector<int>   SubptInd;
	int freqTot = getSubPointIndex(filename1,filename2,SubPointIndex);				///SubptInd contains the index of the subsampled points in the original file
	//getchar()
	int count1=0;
	

	for(int i=0;i<freqTot;i++)
		SubptInd.push_back(0);

	typedef std::map<int,int>::iterator it_type;
for(it_type iterator = SubPointIndex.begin(); iterator != SubPointIndex.end(); iterator++) 
		{
			//if(iterator->second > freqTot)
			//	SubptInd[iterator->first] = 0;	
			SubptInd[iterator->first] = iterator->second;
			count1++;
		}
cout<<"\ninside inputGraph. SubptInd: "<<SubptInd.size(); 	
//	ofstream fileee("SubptInd");

	int max1 = 0, max2 = 0;
/*
	for(int i=0;i<SubptInd.size();i++)
			{
				fileee<<SubptInd[i]<<"\n";
				if(SubptInd[i] > SubptInd.size())
					SubPointIndex[i] = 0;
				else if(SubptInd[i]!=0)
					{	SubPointIndex[i] = SubptInd[i];	count1++;}
				else if(SubptInd[i]==0)
					{SubPointIndex[i] = 0;}
				if(SubptInd[i]>max1)
					max1 = SubptInd[i];
			}
*/			//getchar();
	//exit(0);
	cout<<" SubPointIndex "<<SubPointIndex.size()<<"\n"<<" count "<<count1<<".		";
	file>>str;		//get dimension
	dim = atoi(str);
	
	file>>str;	//get number of points
	noPoints = atoi(str);
	string nouse;
	if(dim!=0)					//point coordinates are given
	{
		for(int i=0;i<dim*noPoints;i++)
			file>>nouse;

	}	//all datapoints have been removed. now form the graph
	
	//////////Read in graph
	float *MortonGraph[noPoints];
	
    for (int i=0; i<noPoints; i++)
        {
        	MortonGraph[i] = (float *)malloc(noPoints * sizeof(float));
        	for(int j=0;j<noPoints;j++)
        		{
        			MortonGraph[i][j] = 0.0;
        			//MortonGraph[j][i] = INT_MAX;
        		}

        } 

        //getchar();	
        //memset(MortonGraph, 1000, sizeof(int) * noPoints * noPoints);

          char str1[10], str2[10], str3[10];
     
     while(!file.eof())
     {
     	file>>str1;	file>>str2;	file>>str3;
     	int u = atoi(str1), v = atoi(str2);
     	float weight = atof(str3);
     	MortonGraph[u][v] = weight;
     	MortonGraph[v][u] = weight;

     }

    //cout<<" entering colormap loop.";
	//getchar();
float first[noPoints];
int count=0;
//getchar();
for(int i=0;i<noPoints;i++)
{
	if(SubptInd[i]!=0)
		{ 
			
			if(SubptInd[i]>freqTot)
				ColorMapping[i]= 0;
			else
				ColorMapping[i]= SubptInd[i];

				count++;
		}	//cout<<"end if part "; getchar();}	//this point is selected
	else
	{
		//cout<<"in else portion";
		//getchar();
		//cout<<" if part "; getchar(); 
		int minpos = dijkstra(MortonGraph,noPoints,i,first,SubptInd);
		//cout<<"out of dijkstra";
		//getchar();
		//int minpos = findmin(first,noPoints,SubptInd);
		if(minpos == -1)	{cout<<"Error ! Point no:"<<i; exit(0);}
		//cout<<"before color map";
		//getchar();
		else if(minpos==freqTot+10)
			ColorMapping[i] = 0;
		else if(minpos>SubptInd.size())
			{ColorMapping[i] = 0;	cout<<"Error. greater than subpt index "<<minpos<<"fir index "<<i; exit(0);}
		else
			ColorMapping[i] = minpos;
		//cout<<"after color map";
		//getchar();
	}

}
//cout<<"		exit colormap loop.";
ofstream filee("ColorMap");
//cout<<"color map"<<count<<"\n";
//getchar();
for(int i=0;i<ColorMapping.size();i++)
		{
			filee<<ColorMapping[i]<<" "<<i<<" \n";
			if(ColorMapping[i]>max2)
				max2 = ColorMapping[i];
			//if(SubptInd[i]==0)

		}

cout<<"color map. count: "<<count<<" max1: "<<max1<<" max2: "<<max2<<"SubPointIndex size from inputGraph: "<<SubPointIndex.size()<<"\n";
	
//cout<<"returned";
//getchar();

//cout<<"returned";
//getchar();





  ////////////////Program over. delete morton code
for (int i=0; i<noPoints; i++) 
        delete MortonGraph[i];

    //delete [] MortonGraph;

}	//end of function inputGraph




//////////////////////////////
void ConstructRipsGraph(const float eps_val, const PointSet &inPts,  SimpleGraph &outRipsGraph)
{
//	Construct Rips graph on the input points
	float sqDistance = eps_val * eps_val;
	ANNSearch::BuildDistanceGraph_ann(outRipsGraph, inPts, sqDistance);
	return;
}
void ComputeEdgeWeights(const PointSet &inPts, SimpleGraph &outRipsGraph)
{
	float segLength = 0.f;
	LongVector tempDiff;
	//
	for (unsigned int i = 0; i < outRipsGraph.vecNode.size(); i++)
	{
		for (std::map<int, float>::iterator mIter = outRipsGraph.vecNode[i].edgeWeights.begin();
				mIter != outRipsGraph.vecNode[i].edgeWeights.end();
				mIter++)
		{
			tempDiff = (*inPts._PointSet)[i] - (*inPts._PointSet)[mIter->first];
			segLength = norm(tempDiff);
			mIter->second = segLength;
		}
	}
	return;
}
void DeltaSparseSampling_EuclideanDistance(const PointSet &inPts,  const float delta_dist,
					std::map<int, int> &SubPointIndex,
					std::vector<int> &ColorMapping)
{
	//
	float sqDistance = 0.0;
	//
	// Perform the delta-sampling delta-sparse subsampling
	//
	sqDistance = delta_dist * delta_dist;
	ColorMapping.resize(inPts._PointSet->size());
	ANNSearch::Subsampling_EuclideanDistance(sqDistance,
											inPts,
											SubPointIndex,
											ColorMapping);

	return;
}
void DeltaSparseSampling_GraphDistance( const SimpleGraph &eps_rips_graph,
										const float delta_dist,
										std::map<int, int> &SubPointIndex,
										std::vector<int> &ColorMapping)
{
	//
	// Perform the delta-sampling delta-sparse subsampling
	//
	ColorMapping.resize(eps_rips_graph.vecNode.size());
	ANNSearch::Subsampling_GraphDistance(delta_dist,
										 eps_rips_graph,
										 SubPointIndex,
										 ColorMapping);
	return;
}
/////////////////

void ComputeBiColoredGraph(const float eps_sampling_dist, const float delta_dist, const bool graph_distance_flag, const PointSet &inPts,
						std::map<int, int> &SubPointIndex, SimpleGraph &biColoredGraph)
{
	//
	SimpleGraph RipsGraph;
	std::vector<int> ColorMapping;
	//
	//getSubPointIndex(SubPointIndex,inPts);
	//cout<<"\n\n";
	ConstructRipsGraph(eps_sampling_dist, inPts, RipsGraph);

	//
	//std::cout << "org comp " << RipsGraph.CheckComponents() << std::endl;
	// generate the batch file for creating rips using this rips graph
	// Delat-Sparse-Sampling and construct rips graph on sub-sampling points
	if (graph_distance_flag)
	{
		ComputeEdgeWeights(inPts, RipsGraph);
		//
//		std::vector<std::vector<float> > pts;
//        		std::vector<float> tmp(inPts._dimension);
//        		pts.reserve(inPts._PointSet->size());
//        		for (unsigned int i = 0; i < inPts._PointSet->size(); i++)
//        		{
//
//        			for (unsigned int j = 0; j < inPts._dimension; j++)
//        			{
//        				tmp[j] = (*inPts._PointSet)[i][j];
//        			}
//        			pts.push_back(tmp);
//        		}
//        		RipsGraph.WriteWeightedGraph("test_wg_pts.txt", pts);
//        		exit(0);
		//
		//RipsGraph.WriteWeightedGraph("test_wg.txt", pts);
		//exit(0);
		//
		DeltaSparseSampling_GraphDistance(RipsGraph, delta_dist, SubPointIndex, ColorMapping);
	}
	else
		DeltaSparseSampling_EuclideanDistance(inPts, delta_dist, SubPointIndex, ColorMapping);
	//
	//for(int i=0;i<ColorMapping.size();i++)
	//	cout<<ColorMapping[i]<<" "<<i<<" \n";
	// construct the bi-colored graph from the org rips graph
	//cout<<"Bi-color sub points\n";



	for(int i =0;i<SubPointIndex.size();i++)
		cout<<SubPointIndex[i]<<"  ";

	

	ANNSearch::SetColorMappingAndExtractColoredGraph(ColorMapping, RipsGraph, biColoredGraph);
	//
	//std::cout << "bic comp " << biColoredGraph.CheckComponents() << std::endl;
	biColoredGraph.color_number = (int)SubPointIndex.size();
	//RipsGraph.WriteBackToFile("test_rip.txt");
	return;
}

void ComputeBiColoredGraph(const char* pFileName, const float delta_dist, std::vector<std::vector<float> > &pts,
						std::map<int, int> &SubPointIndex, SimpleGraph &biColoredGraph,const char* filename2)
{
	//
	//cout<<"till here";			getchar();
	SimpleGraph orgGraph;
	std::vector<int> ColorMapping;
	//
	//cout<<"till here2";			getchar();
	orgGraph.ReadWeightedGaph(pFileName, pts);

	//cout<<"came in";	getchar();

	ColorMapping.resize(orgGraph.vecNode.size());
	inputGraph(pFileName, filename2, SubPointIndex, ColorMapping );
//	exit(0);
	//
	//std::cout << "org comp " << RipsGraph.CheckComponents() << std::endl;
	// generate the batch file for creating rips using this rips graph
	// Delat-Sparse-Sampling and construct rips graph on sub-sampling points

	//DeltaSparseSampling_GraphDistance(orgGraph, delta_dist, SubPointIndex, ColorMapping);
	//std::cout<<" after all this modified. SubPointIndex:"<<SubPointIndex.size()<<" ColorMapping: "<<ColorMapping.size()<<"\n";
	//ofstream fileee("SubptInd2");
	//for(int i =0;i<SubPointIndex.size();i++)
	//	fileee<<SubPointIndex[i]<<" ||||| \n";

	//
	// cout<<"Bi-color\n";
	//for(int i =0;i<ColorMapping.size();i++)
	//	cout<<ColorMapping[i]<<" "<<i<<"\n";
	// construct the bi-colored graph from the org rips graph
	//cout<<"went out";
	//getchar();
	ANNSearch::SetColorMappingAndExtractColoredGraph(ColorMapping, orgGraph, biColoredGraph);
	ofstream filee("ColorMap");
	ofstream fileee("SubPointIndex");
	int max2 = 0;
	int max3 = 0;

typedef std::map<int,int>::iterator it_type;
for(it_type iterator = SubPointIndex.begin(); iterator != SubPointIndex.end(); iterator++) 
{
	//if(iterator->second!=0)
	{
	fileee<< iterator->first<<" ";
    fileee<< iterator->second<<"\n" ;
	}
    
    // Repeat if you also want to iterate through the second map.
}
	//cout<<"color map"<<count<<"\n";
	//getchar();
	for(int i=0;i<ColorMapping.size();i++)
		{
			filee<<ColorMapping[i]<<" "<<i<<"\n";
			if(ColorMapping[i]>max2)
			 	max2 = ColorMapping[i];
			//if(SubptInd[i]==0)

		}

	//
	//std::cout << "bic comp " << biColoredGraph.CheckComponents() << std::endl;
	//std::cout<<" after all this modified 2. SubPointIndex:"<<SubPointIndex.size()<<"\n";
	biColoredGraph.color_number = (int)SubPointIndex.size();
	std::cout<<" after all this modified 3. SubPointIndex:"<<SubPointIndex.size()<<" ColorMapping max:"<<max2<<" ColorMapping size:"<<ColorMapping.size()<<"\n";
	//RipsGraph.WriteBackToFile("test_rip.txt");
	return;
}
void WriteOFFformatComplex(const char* pFileName, const std::vector<std::vector<float> > &pts, std::map<int, int> &subPointIndices, GIComplex &gic)
{
	std::cout << "Writing < " << pFileName << " >" << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		if (!pts.empty())
		{
			// write points out
			std::vector<int> vecSubPts(subPointIndices.size());
			//cout<<"from writing part: SubPointIndex size:"<<subPointIndices.size();
			for (std::map<int, int>::iterator mIter = subPointIndices.begin();
				mIter != subPointIndices.end(); mIter++)
			{
				vecSubPts[mIter->second] = mIter->first;
			}
			//
			sstr << pts.front().size() << " " << vecSubPts.size() << std::endl;
			for (unsigned int i = 0; i < vecSubPts.size(); i++)
			{
				for (unsigned v = 0; v < pts.front().size(); v++)
				{
					sstr << pts[vecSubPts[i]][v] << " ";
				}
				sstr << std::endl;
			}
		}
		else
		{
			sstr << "0 " << subPointIndices.size() << std::endl;
			std::vector<int> subVertexIndices(subPointIndices.size());
			for (std::map<int, int>::iterator mIter = subPointIndices.begin();
				mIter != subPointIndices.end(); mIter++)
				subVertexIndices[mIter->second] = mIter->first;
			for (unsigned int i = 0; i < subVertexIndices.size(); i++)
				sstr << subVertexIndices[i] << std::endl;
		}
		//
		// now add edges
		//
		for (int i = 1; i < gic.dim + 1; i++)
		{
			//depth = i + 1; // want to visit all simplices at this dimension
			if (gic.Simplicies_Cnt[i] != 0)
			{// visit each simplex in dimension i
				for (int vid = 0; vid < gic.Simplicies_Cnt[0]; vid++)
				{
					if (gic.head_circular_list_in_each_dim[vid][i - 1]) // i-1 as vertex are stored in an array
					{// the circular list is not empty
						SimplicialTreeNode_ptr pIter(gic.head_circular_list_in_each_dim[vid][i - 1]);
						do
						{// visit each simplex
							sstr << i + 1<< " ";
							SimplicialTreeNode_ptr pParentIter(pIter);
							//
							do
							{
								sstr << pParentIter->last_label << " ";
								//
								pParentIter = pParentIter->parent_ptr;
							}while (pParentIter);
							sstr << std::endl;
							// move ahead
							pIter = pIter->next_circular_ptr;
						}while (pIter != gic.head_circular_list_in_each_dim[vid][i - 1]);
					} // if head
				}// for vid
			} // if si
		}// for i dim
		sstr << std::endl;
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "---- Done ---- " << std::endl << std::endl;
	//
}
void WriteOFFformatComplex(const char* pFileName, const PointSet &pts, std::map<int, int> &subPointIndices, GIComplex &gic)
{
	std::cout << "Writing < " << pFileName << " >" << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		// write points out
		std::vector<int> vecSubPts(subPointIndices.size());
		//std::cout<<"after all this: SubPointIndex: ";//<<subPointIndices.size()<<"\n";
		for (std::map<int, int>::iterator mIter = subPointIndices.begin();
			mIter != subPointIndices.end(); mIter++)
		{
			vecSubPts[mIter->second] = mIter->first;
		}
		//
		sstr << pts._dimension << " " << vecSubPts.size() << std::endl;
		for (unsigned int i = 0; i < vecSubPts.size(); i++)
		{
			for (unsigned v = 0; v < pts._dimension; v++)
			{
				sstr << (*pts._PointSet)[vecSubPts[i]][v] << " ";
			}
			sstr << std::endl;
		}
		//
		// now add edges
		//
		for (int i = 1; i < gic.dim + 1; i++)
		{
			//depth = i + 1; // want to visit all simplices at this dimension
			if (gic.Simplicies_Cnt[i] != 0)
			{// visit each simplex in dimension i
				for (int vid = 0; vid < gic.Simplicies_Cnt[0]; vid++)
				{
					if (gic.head_circular_list_in_each_dim[vid][i - 1]) // i-1 as vertex are stored in an array
					{// the circular list is not empty
						SimplicialTreeNode_ptr pIter(gic.head_circular_list_in_each_dim[vid][i - 1]);
						do
						{// visit each simplex
							sstr << i + 1<< " ";
							SimplicialTreeNode_ptr pParentIter(pIter);
							//
							do
							{
								sstr << pParentIter->last_label << " ";
								//
								pParentIter = pParentIter->parent_ptr;
							}while (pParentIter);
							sstr << std::endl;
							// move ahead
							pIter = pIter->next_circular_ptr;
						}while (pIter != gic.head_circular_list_in_each_dim[vid][i - 1]);
					} // if head
				}// for vid
			} // if si
		}// for i dim
		sstr << std::endl;
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "---- Done ---- " << std::endl << std::endl;
	//
}
/***********************************************/
char * strLicense = "THIS SOFTWARE IS PROVIDED \"AS-IS\". THERE IS NO WARRANTY OF ANY KIND. "
"NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR "
"ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY. \n"
"\n"
"This software was developed (and is copyrighted by) the Jyamiti group at "
"The Ohio State University. Please do not redistribute this software. "
"This program is for academic research use only. This software uses the "
"Boost library (www.boost.org) and Ann library "
"(www.cs.umd.edu/~mount/ANN/) which are covered under their own licenses.\n"
"\n"
"The Boost library's license "
"(which applies to the Boost library ONLY and NOT to this program itself) is "
"as follows:\n"
"\n"
"LICENSE\n"
"---------------------------------------------------------------------------\n"
"Boost Software License - Version 1.0 - August 17th, 2003\n"
"\n"
"Permission is hereby granted, free of charge, to any person or organization "
"obtaining a copy of the software and accompanying documentation covered by "
"this license (the \"Software\") to use, reproduce, display, distribute, "
"execute, and transmit the Software, and to prepare derivative works of the "
"Software, and to permit third-parties to whom the Software is furnished to "
"do so, all subject to the following: \n"
"\n"
"The copyright notices in the Software and this entire statement, including "
"the above license grant, this restriction and the following disclaimer, "
"must be included in all copies of the Software, in whole or in part, and "
"all derivative works of the Software, unless such copies or derivative "
"works are solely in the form of machine-executable object code generated by "
"a source language processor. \n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "
"FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT "
"SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE "
"FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, "
"ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER "
"DEALINGS IN THE SOFTWARE. \n"
"---------------------------------------------------------------------------\n"
"\n"
"ANN library's license "
"(which applies to the ANN library ONLY and NOT to this program itself) is "
"as follows: \n"
"\n"
"LICENSE\n"
"---------------------------------------------------------------------------\n"
"The ANN Library (all versions) is provided under the terms and "
"conditions of the GNU Lesser General Public Library, which is stated "
"below.  It can also be found at: \n"
"\n"
"   http:////www.gnu.org/copyleft/lesser.html \n"
"---------------------------------------------------------------------------\n";
/**********************************************************************/
bool ParseCommand(int argc, char** argv,
				int &input_type,
				std::string &InputFile,
				std::string &SubsampleFile,
				std::string &OutputFile,
				bool &graph_distance_flag,
				float &eps_dist,
				float &delta_dist,
				int &max_dimension)
{
	try
	{
		/* Define the program options description
		*/
		namespace po = boost::program_options;
		po::options_description desc("GIComplex Usage");
		desc.add_options()
			(",h", "Help information;")
			(",l", "License information;")
			("type", po::value<int>(&input_type)->required(), "Input data type: point cloud (0) or weighted graph (1);")
			(",I", po::value<std::string>(&InputFile)->required(), "Input file name;")
			(",S", po::value<std::string>(&SubsampleFile)->required(), "Subsampled file name;")
			(",O", po::value<std::string>(&OutputFile)->required(), "Output file name prefix;")
			(",g", po::value<bool>(&graph_distance_flag)->default_value(false), "Subsample by Euclidean distance (false) or graph distance (true); for weighted graph input, it's always (true);")
			("epsilon", po::value<float>(&eps_dist), "Parameter for building the Rips graph on input samples; Not used for weighted graph input;")
			("delta", po::value<float>(&delta_dist)->required(), "Parameter for generating subsamples from input samples;")
			("dim", po::value<int>(&max_dimension)->required(), "The maximum dimension of simplices in output complex;");

		// Parser map
		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, desc), vm);

			//
			if (vm.count("-h"))
			{
				std::cout << desc << std::endl;
			}
			//
			if (vm.count("-l"))
			{
				std::cout << strLicense << std::endl;
			}
			//
			po::notify(vm);
		}
		catch(boost::program_options::required_option& e)
		{
			std::cerr<< "ERROR: " << e.what() << std::endl;
			return false;
		}
		catch(boost::program_options::error& e)
		{
			std::cerr<< "ERROR: " << e.what() << std::endl;
			return false;
		}
	}
	catch(std::exception& e)
	{
		std::cerr << "Unhandled Exception reached the top of main: "
					<< e.what() << ", application will now exit" << std::endl;
		return false;

	}
	return true;
}
/*********************/
int main(int argc, char **argv)
{
	//
	std::string InputFileName;
	std::string OutputFileName;
	std::string SubsampleFileName;
	int input_type = 0;
	float eps_dist = .05f;
	float delta_dist = 0.9f;
	bool graph_distance_flag = true;
	int max_dimension = 3;
	std::cout<<"Beware ! Hacked version !!\n \n || Remember to change format of subsample file? ||\n\n || Remember to change Bill's graph code to Gic format? ||\n";
	//getchar();
//	inputGraph("graph.txt","datapoints-subsample.txt");
//	cout<<"\nfinish read";
//	exit(0);
	//
	//std::cout<<"graphBunny.txt ";
	if (ParseCommand(argc, argv, input_type, InputFileName, SubsampleFileName, OutputFileName, graph_distance_flag, eps_dist, delta_dist, max_dimension))
	{
		if (!input_type)
		{
			PointSet pts;

			//pts.SampleCrossLine(0.02, 1.0);
			pts.ReadPointsFromFile(InputFileName.c_str());
			std::map<int, int> SubPointIndex;
			SimpleGraph biColoredGraph;
			//
			std::cout << std::endl << "Construct graph induced complex " << std::endl << std::endl;
			//
			ComputeBiColoredGraph(eps_dist, delta_dist, graph_distance_flag, pts, SubPointIndex, biColoredGraph);
			//
			
			std::string outfile_name(OutputFileName);
			//
			GIComplex gic(max_dimension, biColoredGraph.color_number, &biColoredGraph);
			gic.Construction();
			//
			
			outfile_name = outfile_name + "_stat.txt";
			gic.WriteStatisticsToFile(outfile_name.c_str());
			//
			outfile_name = OutputFileName + "_complex.txt";
			WriteOFFformatComplex(outfile_name.c_str(), pts, SubPointIndex, gic);
		}
		else
		{// read weighted graph
			std::vector<std::vector<float> > pts;
			std::map<int, int> SubPointIndex;
			SimpleGraph biColoredGraph;
			//

			ComputeBiColoredGraph(InputFileName.c_str(), delta_dist, pts, SubPointIndex, biColoredGraph,SubsampleFileName.c_str());
			//std::cout<<"from main SubPointIndex "<<biColoredGraph.color_number;
			//
			//
			std::cout << std::endl << "Construct graph induced complex " << std::endl << std::endl;
			//
			std::string outfile_name(OutputFileName);
			//
			GIComplex gic(max_dimension, biColoredGraph.color_number, &biColoredGraph);
			//cout<<"here";
			//getchar();
			gic.Construction();
			//
			//getchar();
			std::cout<<"from main SubPointIndex : "<<SubPointIndex.size();
			outfile_name = outfile_name + "_stat.txt";
			gic.WriteStatisticsToFile(outfile_name.c_str());
			//
			outfile_name = OutputFileName + "_complex.txt";
			//getchar();
			WriteOFFformatComplex(outfile_name.c_str(), pts, SubPointIndex, gic);
		}
		//pts.WriteBackToFile("test.txt");
		//biColoredGraph.WriteBackToFile("test_graph.txt");
	}
	//
	return 0;
}
