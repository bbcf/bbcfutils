#include <stdio.h>
#include <cstring>
#include <string>
#include <math.h>
//#include <glib.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include "slidetools.h"
#include "FixSizeSortedList.h"
// #define CUTOFF 0.0f
// #define CUTOFF -1.0e12
#define MAX_NAME 256

using namespace std;

char* reverse(char s[])
{
  int c, i, j;
   
  for (i=0, j=strlen(s)-1; i < j;++i, --j)
  {
    c = s[i];
    s[i] = s [j];
    s[j] = c;
  }
  return s;
}

string* reverse(string* s)
{
  int c, i, j;
   
  for (i=0, j=s->size()-1; i < j;++i, --j)
  {
    c = (*s)[i];
    (*s)[i] = (*s)[j];
    (*s)[j] = c;
  }
  return s;
}

string* complement(string* s)
{
  for (int i = 0; i < s->size();++i){
    switch((*s)[i]){
    case 'a':
      {(*s)[i]='t';break;}
    case 'A':
      {(*s)[i]='T';break;}
    case 'c':
      {(*s)[i]='g';break;}
    case 'C':
      {(*s)[i]='G';break;}
    case 'g':
      {(*s)[i]='c';break;}
    case 'G':
      {(*s)[i]='C';break;}
    case 't':
      {(*s)[i]='a';break;}
    case 'T':
      {(*s)[i]='A';break;}
    } 
  }
  return s;
}

void addDashs(string* s, int n)
{
  while(n>0){
    *s+="_";
    --n;
  }
}

int main(int argc, char **argv)
{
  double cutoff;

  char buf[BUF];
  FILE *fin=NULL;

  WeightMat *wm1=NULL, *bg=NULL, *wm1R=NULL, *wm2=NULL, *wm2R=NULL;
  cerr << "Number of arguments = " << argc << endl;
  if(argc<7)
  {
    fprintf(stderr, "usage: S2M wm1 wm2 bg spacer cutoff fasta\n");
    exit(0);
  }

  wm1=read_WM(argv[1]); //Read weight matrix
  wm1R=read_WMR(argv[1]); //Read weight matrix
  wm2=read_WM(argv[2]); //Read weight matrix
  wm2R=read_WMR(argv[2]); //Read weight matrix
  bg=read_WM(argv[3]); //Read background matrix
  const int maxMatSize=max(wm1->span,wm2->span);
  const int gap=max(1,atoi(argv[4])+1);
  cerr << "GAP = " << gap-1 << endl;
  cerr << "Size Matrix 1 = " << wm1->span << endl;
  cerr << "Size Matrix 2 = " << wm2->span << endl;

  fsslist bestPrevious(gap);
  fsslist bestPreviousR(gap);
  list<seqScore> buffer;
  list<seqScore> bufferR;

  {
	istringstream is(argv[5]);
	is >> cutoff;
  }

  cerr << "Score cutoff = " << cutoff << endl;

  fin=fopen(argv[6], "r");
  cerr << "Filename = " << argv[6] << endl;

  int readed=0;
  int state=0; //0=Waiting >,1=Parse seq name,2=Scan W Matrix
  int renduName=0;
  int positionInSequence=0;
  char name[MAX_NAME+1];
  int offset=0;
  bool reachSpace=false;

  while ((readed=fread(buf+offset,sizeof(char),BUF-offset,fin))){
    readed+=offset;
    offset=0;
    for (int currentPos=0;currentPos < readed;++currentPos){
      switch(state){
      case 0:{
	if (buf[currentPos] == '>'){
	  state=1;
	  renduName=0;
	}
	break;
      }
      case 1:{
	if (buf[currentPos] == '\n'){
 
	  name[renduName]='\0';
	  cerr << name << endl;
	  //printf("%s\n",name);
	  state=2;
	  renduName=0;
	  positionInSequence=1;
	  reachSpace=false;
	  bestPrevious.clear();
	  bestPreviousR.clear();
	  buffer.clear();
	  bufferR.clear();
	}
	else if (renduName < MAX_NAME && !reachSpace){
	  if (buf[currentPos] != ' '){
	    name[renduName++]=buf[currentPos];
	  }
	  else{
	    reachSpace=true;
	  }
	}
	break;
      }
      case 2:{
	if (buf[currentPos] == '>'){
	  state=1;
	}
	else if (buf[currentPos] != '\n'){
	  //cerr << positionInSequence  << endl;
	  double score1=0.0;
	  double score1R=0.0;
	  int posinMat1=0;
	  char curSeq1[wm1->span+1];curSeq1[0]='\0';
	  int renduSeq1=0;
	  for(renduSeq1=0; posinMat1 < wm1->span && currentPos+renduSeq1 < readed; ++renduSeq1){
	    if(wm1->pos[posinMat1].flag){
	      if (buf[currentPos+renduSeq1] == '>'){
		state=0;
		break;
	      }
	      else if(buf[currentPos+renduSeq1] == '\n'){//Skip this char but the matrix stay at the same place
	      }
	      else{
		curSeq1[posinMat1]=buf[currentPos+renduSeq1];
		score1 += wm1->pos[posinMat1].f[buf[currentPos+renduSeq1]] - 
		  bg->pos[0].f[buf[currentPos+renduSeq1]];
		score1R += wm1R->pos[posinMat1].f[buf[currentPos+renduSeq1]] - 
		  bg->pos[0].f[buf[currentPos+renduSeq1]];
		
		++posinMat1;
	      }
	    }
	  }

	  double score2=0.0;
	  double score2R=0.0;
	  int posinMat2=0;
	  char curSeq2[wm2->span+1];curSeq2[0]='\0';
	  int renduSeq2=0;
	  for(renduSeq2=0; posinMat2 < wm2->span && currentPos+renduSeq2 < readed; ++renduSeq2){
	    if(wm2->pos[posinMat2].flag){
	      if (buf[currentPos+renduSeq2] == '>'){
		state=0;
		break;
	      }
	      else if(buf[currentPos+renduSeq2] == '\n'){//Skip this char but the matrix stay at the same place
	      }
	      else{
		curSeq2[posinMat2]=buf[currentPos+renduSeq2];
		score2 += wm2->pos[posinMat2].f[buf[currentPos+renduSeq2]] - 
		  bg->pos[0].f[buf[currentPos+renduSeq2]];
		score2R += wm2R->pos[posinMat2].f[buf[currentPos+renduSeq2]] - 
		  bg->pos[0].f[buf[currentPos+renduSeq2]];
		
		++posinMat2;
	      }
	    }
	  }

	  bool doesntReachTheEndForTheLongestMatrix = false;
	  if (posinMat1 < wm1->span &&
	      wm1->span == maxMatSize){
	    doesntReachTheEndForTheLongestMatrix=true;
	  }
	  else if(posinMat2 < wm2->span &&
		  wm2->span == maxMatSize){
	    doesntReachTheEndForTheLongestMatrix=true;
	  }

	  if (doesntReachTheEndForTheLongestMatrix){
	    if(state == 2){ //we don't reach the end of the matrix and
	      //we are still in state 2!
	      for(offset=0; currentPos+offset < readed; ++offset){
		buf[offset]=buf[currentPos+offset];
	      }
	      currentPos=readed;//Read another buffer
	    }
	  }
	  else{
	    curSeq1[posinMat1]='\0';
	    curSeq2[posinMat2]='\0';
	    //cerr << curSeq1 << "\t" << curSeq2 << endl;

	    //positionInSequence == 1 ?
	    //  bestPrevious.insert(make_pair(score1,make_pair(positionInSequence,new string(curSeq1)))) :
	    buffer.push_back(make_pair(score1,make_pair(positionInSequence,new string(curSeq1))));

	    //positionInSequence == 1 ?
	    //  bestPreviousR.insert(make_pair(score2R,make_pair(positionInSequence,new string(curSeq2)))) :
	    bufferR.push_back(make_pair(score2R,make_pair(positionInSequence,new string(curSeq2))));

	    if(positionInSequence>wm1->span){
	      bestPrevious.insert(buffer.front());
	      buffer.pop_front();
	      seqScore maxPrevious=bestPrevious.max();
	      
	      if (maxPrevious.first+score2 > cutoff){
		string temp=*(maxPrevious.second.second);
		addDashs(&temp,positionInSequence-(maxPrevious.second.first+wm1->span));
		temp+=curSeq2;
		cout << name << "\t"
		     << temp << "\t"
		     << maxPrevious.first+score2 << "\t"
		     << maxPrevious.second.first << "\t+" << endl;
	      }
	    }

	    if(positionInSequence>wm2->span){
	      bestPreviousR.insert(bufferR.front());
	      bufferR.pop_front();
	      seqScore maxPrevious=bestPreviousR.max();
	      if (maxPrevious.first+score1R > cutoff){
		string temp=*(maxPrevious.second.second);
		addDashs(&temp,positionInSequence-(maxPrevious.second.first+wm2->span));
		temp+=curSeq1;
		reverse(complement(&temp));
		cout << name << "\t"
		     << temp << "\t"
		     << maxPrevious.first+score1R << "\t"
		     << maxPrevious.second.first << "\t-" << endl;
	      }
	    }
	    ++positionInSequence;
	  }
	}
	break;
      }
      }
    }
  }
  return 0;
}
