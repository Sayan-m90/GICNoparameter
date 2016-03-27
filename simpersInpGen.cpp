#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>	
#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <ANN/ANN.h>
#include "sampler.h"

using namespace std;
int timeee;
int vertices=0, edges=0 , faces=0 , quads=0;


void Sortfour(int &a, int &b, int &c, int &d)
{
	int low1,high1,low2,high2,middle1,middle2,highest,lowest;

    if ( a < b )
    {
    	low1 = a;
        high1 = b;
    }
        
    else 
        {
        	low1 = b;
        	high1 = a;
        }

    if (c < d)
        {
        	low2 = c;
        	high2 = d;
    	}	
    else
        {
        	low2 = d;
        	high2 = c;
        }

    if (low1 < low2)
        {
        	lowest = low1;
        	middle1 = low2;
        }
    else
        {
        	lowest = low2;
        	middle1 = low1;
        }

    if (high1 > high2)
        {
        	highest = high1;
        	middle2 = high2;
        }
    else

        {
        	highest = high2;
        	middle2 = high1;
        }

    if (middle1 < middle2)
    {
    	a = lowest;
    	b = middle1;
    	c = middle2;
    	d = highest;
        //return (lowest,middle1,middle2,highest);
    }
    else
    {
    	a = lowest;
    	b = middle2;
    	c = middle1;
    	d = highest;
        //return (lowest,middle2,middle1,highest);
    }

}


void Sort(int &a, int &b, int &c){
    if(a>b){
        int tmp = a;
        a = b;
        b = tmp;
    }
    if(a>c){
        int tmp = a;
        a=c;
        c = tmp;
    }
    if(b>c){
        int tmp = b;
        b=c;
        c=tmp;
    }
    return;
}

void append(FILE *head, FILE *tail)
{
    char buff[50];
    size_t n;
    int count = 0;
    while(fscanf(tail,"%s",buff)!=EOF)
    {
    	

    	if(strcmp(buff,"||")==0)
				fprintf(head,"# %d\n",++timeee*10);
		else
				{
					fprintf(head, "%s ", buff);
					count++;
				}
	if(count == 4)
		{
			fprintf(head, "\n");		
			count = 0;
		}
    }

}

int vertexcount()
{
	FILE *fp= fopen("outGIC_stat.txt","r");
	char buff[100];
	int i;

	for(i=0;i<3;i++)
		{
			fscanf(fp,"%s",buff);
			cout<<buff<<" ";
		}

	int vert;
	fscanf(fp,"%d",&vert);
	fclose(fp);
	return vert;

}

int callappend()
{
	FILE	*file1 = fopen("inpSimp.txt","ab");
FILE	*file2 = fopen("collapsecoord.txt","rb");

append(file1,file2);
//fprintf(file1,"# %d\n",++timeee*0);
		//timeee++;
//fprintf(file1,"#5\n");

fclose(file1);
fclose(file2);		

}

void checkstat()
{
	FILE *fs = fopen("outGIC_stat.txt","r");
	char buff[50];
	fscanf(fs,"%s",buff);		// DIM
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// 0
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// :
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// #vertices
	vertices = atoi(buff);
	cout<<vertices<<" \n";

	fscanf(fs,"%s",buff);		// DIM
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// 1
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// :
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// #vertices
	edges = atoi(buff);
	cout<<edges<<" \n";

	fscanf(fs,"%s",buff);		// DIM
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// 2
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// :
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// #vertices
	faces = atoi(buff);
	cout<<faces<<" \n";

	fscanf(fs,"%s",buff);		// DIM
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// 3
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// :
	cout<<buff<<" ";
	fscanf(fs,"%s",buff);		// #vertices
	quads = atoi(buff);
	cout<<quads<<" \n";
	fclose(fs);
}
int main(int argc, char **argv)
{
	
	FILE *fp = fopen("outGIC_complex.txt","r");
	FILE *fo = fopen("inpSimp.txt","w");
	
	checkstat();
	
	float **coordinates;			// 3D coordinate of each point
	//coordinates = (float**)malloc(n*sizeof(float*));
	char buff[1000];
	fscanf(fp,"%s",buff);		// dimension of points, useless
	fscanf(fp,"%s",buff);		// no of points
	int n = atoi(buff);
	int i, j, dim, temp;
	int v1,v2;
	int bluff;
	int e1,e2,e3;
	int q1,q2,q3,q4;
	int vert = vertexcount();
	int dim1 = vert, dim2 = vert, dim3 = vert;
	cout<<"vert:"<<vert;
	//getchar();
	int **traing = (int **)malloc(vert*sizeof(int*));

        for (i = 0; i< vert; i++) 
        	traing[i] = (int *) malloc(vert*sizeof(int));

    
    /////////////////////////////////////////////////////////////////////
	///////////////////////////  extracts corrdinates from output of GIC and puts it in a file 'fout'
	////////////////////////// 	Vertex coordinates to be ignored
	////////////////////////////////////////////////////////////////////
	for(i=0;i<n;i++)		
		{
			
			fscanf(fp,"%s",buff);
			
			fscanf(fp,"%s",buff);
			
			fscanf(fp,"%s",buff);
			
			//fprintf(fo, "%s\n", );
		
		}											//3D coordinates scanned
		timeee = 0;
		//fprintf(fo,"# %d\n",++timeee*10);

		////////////// 	insert edges	//////////////////////

		for(i=0;i<edges;i++)		
		{
			fscanf(fp,"%d",&bluff);
			if(bluff != 2)
			{
				cout<<"wrong prgram 2 "<<bluff;
				exit(0);
			}
			
			fscanf(fp,"%d",&v1);
			
			fscanf(fp,"%d",&v2);

			if(v1<v2)
			fprintf(fo,"insert %d %d\n",v1,v2);
			else
			fprintf(fo,"insert %d %d\n",v2,v1);
		
		}
		fprintf(fo,"# %d\n",++timeee*10);

		////////////// 	insert faces	//////////////////////
		for ( i = 0; i < faces; i++)
		{
			fscanf(fp,"%d",&bluff);
			fscanf(fp,"%d",&e1);
			fscanf(fp,"%d",&e2);
			fscanf(fp,"%d",&e3);

			//cout<<"\n"<<bluff<<" "<<e1<<" "<<e2<<" "<<e3<<"\n";
			//getchar();
			
			if(bluff != 3)
			{
				cout<<"\n wrong program 3 "<<bluff;
				exit(0);
			}

			Sort(e1,e2,e3);
			fprintf(fo,"insert %d %d %d\n",e1,e2,e3);			
			
		}

		fprintf(fo,"# %d\n",++timeee*10);

		////////////// 	insert quads	//////////////////////
		for ( i = 0; i < quads; i++)
		{
			fscanf(fp,"%d",&bluff);
			fscanf(fp,"%d",&q1);
			fscanf(fp,"%d",&q2);
			fscanf(fp,"%d",&q3);
			fscanf(fp,"%d",&q4);
			//cout<<"\n"<<bluff<<" "<<e1<<" "<<e2<<" "<<e3<<"\n";
			//getchar();
			
			if(bluff != 4)
			{
				cout<<"\n wrong program 3 "<<bluff;
				exit(0);
			}

			Sortfour(q1,q2,q3,q4);
			fprintf(fo,"insert %d %d %d %d\n",q1,q2,q3,q4);			
			
		}

		fprintf(fo,"# %d\n",++timeee*10);
					fclose(fp);
					fclose(fo);

					callappend();
					fprintf(fo,"# %d\n",++timeee*10);
					//timeee++;
					

}