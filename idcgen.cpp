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
#include <iostream>
using namespace std;
int main()
{
	FILE *fp= fopen("gICgraph1hr_stat.txt","r");
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
	cout<<" "<<vert<<"\n";
	FILE *fout= fopen("iDC.txt","w");
	fprintf(fout,"%d\n",vert);
	for(i=0;i<vert;i++)
	{
		fprintf(fout,"%d\n",i);
	}
	
	return 0;
}