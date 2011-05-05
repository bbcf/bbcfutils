// hit_annotator: scanning sequences with PWM
// Copyright (C) 2011, Eric R Paquet, Bernhard Sonderegger, Felix Naef, EPFL
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*#include "rand.h"*/

#define BUF 1048576

#define WM_SIZE 256
#define MAX_HEADER 256
typedef struct
{
	int flag;
	double f[WM_SIZE];
}
onePos;

typedef struct
{
	int npos;
	int span;
        char *header;
	onePos *pos;
}
WeightMat;

WeightMat *alloc_WM(int n)
{
	WeightMat *wm = (WeightMat*) malloc(sizeof(WeightMat));
	wm->span=n;
	wm->npos=0;
	wm->pos = (onePos*) malloc(n*sizeof(onePos));

	//Initialize all the wm to log2(.25)
	for (int j = 0;j<n;++j){
		for (int i = 0; i < WM_SIZE;++i){
			wm->pos[j].f[i]=log2(.25);
		}
	}

	return wm;
}

WeightMat* read_WM(char *fname)
{
	char buf[BUF], *pt;
	WeightMat *wm=NULL;
	FILE *fin=fopen(fname, "r");
	onePos *position=NULL;
        char header[MAX_HEADER];
	int nbase=0;
	while(fgets(buf, BUF, fin))
	{
	    if (buf[0] != '>') nbase++;
	}
	fclose(fin);

	wm=alloc_WM(nbase);
	fin=fopen(fname, "r");

	position=wm->pos;
	while(fgets(buf, BUF, fin))
	{
	    if (buf[0] == '>'){
		int index = 1;
		
		while (index < MAX_HEADER && buf[index] != '\n'){
		    header[index] = buf[index];
		    index++;
		}
		header[index] = '\0';
		wm->header = header;
	    }
	    else{
		int k, np=0;
		char *pt=buf;
		// printf ("%s\n",buf);
		double freq, min_log2_freq;
		sscanf(pt, "%d %n", &(position->flag), &np);pt+=np;
		
		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		position->f['A'] = position->f['a'] = min_log2_freq = log2(freq);

		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		position->f['C'] = position->f['c'] = log2(freq);
		if (min_log2_freq > position->f['C'])
			min_log2_freq = position->f['C'];

		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		position->f['G'] = position->f['g'] = log2(freq);
		if (min_log2_freq > position->f['G'])
			min_log2_freq = position->f['G'];
			
		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		position->f['T'] = position->f['t'] = position->f['U'] = position->f['u'] = log2(freq);
		if (min_log2_freq > position->f['T'])
			min_log2_freq = position->f['T'];
			
		position->f['N'] = position->f['n'] = position->f['X'] = position->f['x'] = min_log2_freq;
	 	//printf("Forwards: %f\t%f\t%f\t%f\t%f\t%f\n",position->f['A'],position->f['C'],position->f['G'],position->f['T'],position->f['N'],position->f['X']);

		if(position->flag) wm->npos++;
		position++;
	    }
	}
	fclose(fin);
	return wm;
}

WeightMat* read_WMR(char *fname)
{
	char buf[BUF], *pt;
	WeightMat *wm=NULL;
	FILE *fin=fopen(fname, "r");
	onePos *position=NULL;

	int nbase=0;
	while(fgets(buf, BUF, fin))
	{
		nbase++;
	}
	fclose(fin);

	wm=alloc_WM(nbase);
	fin=fopen(fname, "r");

	int matPos=nbase-1;
	while(fgets(buf, BUF, fin))
	{
		int k, np=0;
		char *pt=buf;
		//printf ("%s\n",buf);
		double freq, min_log2_freq;
		sscanf(pt, "%d %n", &(wm->pos[matPos].flag), &np);pt+=np;

		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		wm->pos[matPos].f['T'] = wm->pos[matPos].f['t'] = min_log2_freq = log2(freq);

		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		wm->pos[matPos].f['G'] = wm->pos[matPos].f['g'] = log2(freq);
		if (min_log2_freq > wm->pos[matPos].f['G'])
			min_log2_freq = wm->pos[matPos].f['G'];
		
		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		wm->pos[matPos].f['C'] = wm->pos[matPos].f['c'] = log2(freq);
		if (min_log2_freq > wm->pos[matPos].f['C'])
			min_log2_freq = wm->pos[matPos].f['C'];
		
		sscanf(pt, "%lf %n", &freq, &np);pt+=np;
		wm->pos[matPos].f['A'] = wm->pos[matPos].f['a'] = log2(freq);
		if (min_log2_freq > wm->pos[matPos].f['A'])
			min_log2_freq = wm->pos[matPos].f['A'];
		
		wm->pos[matPos].f['n'] = wm->pos[matPos].f['N'] = wm->pos[matPos].f['x'] = wm->pos[matPos].f['X'] = min_log2_freq;
		//printf("Reverse:  %f\t%f\t%f\t%f\t%f\t%f\n",wm->pos[matPos].f['A'], wm->pos[matPos].f['C'], wm->pos[matPos].f['G'],wm->pos[matPos].f['T'],wm->pos[matPos].f['N'],wm->pos[matPos].f['X']);

		if(wm->pos[matPos].flag) wm->npos++;
		--matPos;
	}
	fclose(fin);

	return wm;
}
