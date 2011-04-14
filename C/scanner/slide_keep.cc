#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sstream>
//#include <glib.h>

#include "slidetools.h"
//#define CUTOFF 10.0f
//#define CUTOFF -1.0e12
//#define CUTOFF 0.0f
#define MAX_NAME 256

char* reverse(char s[]){
    int c, i, j;
     
    for (i=0, j=strlen(s)-1; i < j;++i, --j){
        c = s[i];
        s[i] = s [j];
        s[j] = c;
    }
    return s;
}

char* complement(char s[]){
    for (int i = 0; i < strlen(s);++i){
        switch(s[i]){
        case 'a':
            {s[i]='t';break;}
        case 'A':
            {s[i]='T';break;}
        case 'c':
            {s[i]='g';break;}
        case 'C':
            {s[i]='G';break;}
        case 'g':
            {s[i]='c';break;}
        case 'G':
            {s[i]='C';break;}
        case 't':
            {s[i]='a';break;}
        case 'T':
            {s[i]='A';break;}
        } 
    }
    return s;
}

int main(int argc, char **argv){
    int pos, width;
    double cutoff;

    char buf[BUF];
    FILE *fin=NULL;

    WeightMat *wm=NULL, *bg=NULL, *wmR=NULL;
    if(argc<4){
        fprintf(stderr, "usage: S1K wm1 bg cutoff fasta\n");
        exit(0);
    }

    wm=read_WM(argv[1]);    //Read weight matrix
    wmR=read_WMR(argv[1]);  //Read weight matrix
    bg=read_WM(argv[2]);    //Read background matrix
    {
        std::istringstream is(argv[3]);
        is >> cutoff;
    }
    fin=fopen(argv[4], "r");

    int readed              = 0;
    int state               = 0; //0=Waiting >,1=Parse seq name,2=Scan W Matrix
    int renduName           = 0;
    int positionInSequence  = 0;
    int offset              = 0;
    bool reachSpace         = false;
    char name[MAX_NAME+1];

    while (readed=fread(buf+offset,sizeof(char),BUF-offset,fin)){
        readed  += offset;
        offset  = 0;
        for (int currentPos=0;currentPos < readed;++currentPos){
            switch(state){
                case 0:{
                    if (buf[currentPos] == '>'){
                        state       = 1;
                        renduName   = 0;
                    }
                    break;
                }
                case 1:{
                    if (buf[currentPos] == '\n'){
                        name[renduName]     = '\0'; //printf("%s\n",name);
                        state               = 2;
                        renduName           = 0;
                        positionInSequence  = 1;
                        reachSpace          = false;
                    }
                    else if (renduName < MAX_NAME && !reachSpace){
                        if (buf[currentPos] != ' ')
                            name[renduName++]=buf[currentPos];
                        else
                            reachSpace=true;
                    }
                    break;
                }
                case 2:{
                    if (buf[currentPos] == '>'){
                        state=1;
                    }
                    else if (buf[currentPos] != '\n'){

                        double  score   = 0.0;
                        double  scoreR  = 0.0;
                        int     posinMat= 0;
                        int     renduSeq= 0;
                        char    yo      = 0;
                        char    curSeq[wm->span+1];curSeq[0]='\0';
                        for(renduSeq=0; posinMat < wm->span && currentPos+renduSeq < readed; ++renduSeq){
                            if(wm->pos[posinMat].flag){
                                if (buf[currentPos+renduSeq] == '>'){
                                    state=0;
                                    break;
                                }
                                else if(buf[currentPos+renduSeq] == '\n'){//Skip this char but the matrix stay at the same place
                                    yo=posinMat;
                                }
                                else{
                                    curSeq[posinMat]=buf[currentPos+renduSeq];
                                    score   += wm->pos[posinMat].f[buf[currentPos+renduSeq]]    -   bg->pos[0].f[buf[currentPos+renduSeq]];
                                    scoreR  += wmR->pos[posinMat].f[buf[currentPos+renduSeq]]   -   bg->pos[0].f[buf[currentPos+renduSeq]];
                                    ++posinMat;
                                }
                            }
                        }
                        if (posinMat == wm->span){ //we were able to process a complete chunk!
                            curSeq[posinMat]='\0';
                            if (score > cutoff)
                                printf("%s\t%s\t%f\t%d\t%c\n",name,curSeq,score,positionInSequence,'+');
                            if (scoreR > cutoff){
                                printf("%s\t%s\t%f\t%d\t%c\n",name,complement(reverse(curSeq)),scoreR,positionInSequence,'-');
                            }
                            ++positionInSequence;
                        }
                        else{
                            if(state == 2){ //we don't reach the end of the matrix and
                                //we are still in state 2!
                                for(offset=0; currentPos+offset < readed; ++offset){
                                    buf[offset]=buf[currentPos+offset];
                                }
                                currentPos=readed;//Read another buffer
                            }
                        }
                    }
                    break;
                }
            }
        }
    }
    return 0;
}
