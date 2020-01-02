/* Copyright (C) 2012 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */


/*  Last edited: Jun 23 19:33 1995 (gos) */
/*  date unknown (LaDeana Hillier) */
/*  3/23/04 (Angie Hinrichs): use obs/exp from G-G & F '87*/
/* 
  cpg.c: CpG island finder. Larsen et al Genomics 13 1992 p1095-1107
  "usually defined as >200bp with %GC > 50% and obs/exp CpG >
  0.6". Here use running sum rather than window to score: if not CpG
  at postion i, then decrement runSum counter, but if CpG then runSum
  += CPGSCORE.     Spans > threshold are searched
  for recursively.
  
  filename: fasta format concatenated sequence file: 
*/

/*********************************************************


put in options: -p debug to switch on commented out print statements
                -l length threshold
                -t score threshold
                -s to specify score
		-d score (-ve) if not found : default = -1.
		
Really want it to take pattern file:
     Pattern  score_if_found  score_if_not_found

How best to give #pattern hits, %GC in span, o/e ?
Cleanest:
Provide using separate scripts called using output.

Add -g or -m option for "max gap" as in Qpatch.c - score not increased
if reaches this threshold.  This means that patches which are separated
by more than the threshold are never joined.

*********************************************************/

#include "stdio.h"
#include "math.h"
#include "readseq.h"
#include "string.h"
#include <stdlib.h>

/********************* CONTROLS ********************/
#define CPGSCORE  17   /* so that can compare with old cpg - this
			  had CPGSCORE 27, but allowed score to reach
			  0 once without reporting */
#define MALLOCBLOCK 200000
#define MAXNAMELEN 512  /*"length 0, title > MAXNAMELEN or starts with bad character." if too small */

/*------------------------------------------------------*/
int readSequence (FILE *fil, int *conv, char **seq, char **id, char **desc, int *length) ; /*读序列函数*/

void findspans (FILE *fp, int start, int end, char *seq, char *seqname ) ;

void getstats (int start, int end, char *seq, char *seqname, int *ncpg, int *ngpc, int *ng, int *nc, int *nacga, int *nacgc, int *nacgg, int *nacgt, int *nccga, int *nccgc, int *nccgg, int *nccgt, int *ngcga,int *ngcgc,int *ngcgg,int *ngcgt, int *ntcga,int *ntcgc,int *ntcgg, int *ntcgt) ;
char *join1(char *a, char *b);

void usage (void)
{ 
  fprintf(stderr, "cpglh - calculate CpG Island data for cpgIslandExt tracks\n");
  fprintf(stderr, "usage:\n    cpglh <sequence.fa>\n") ;
  fprintf(stderr, "where <sequence.fa> is fasta sequence, must be more than\n");
  fprintf(stderr, "   200 bases of legitimate sequence, not all N's\n");
  exit (-1);
}
/*------------------------------------------------------*/
int main (int argc, char **argv) /*  主函数 */
{ 
  
  char *seq, *seqname, *desc ;
  
  int length ;
  int i ;
  static FILE *fil ;

  /*------------------------------------------------------*/  
  switch ( argc )
    {
    default: if (argc != 2)
      usage () ;
    }
     if (!(seqname = malloc (MAXNAMELEN+1)))
     { fprintf (stderr, "Couldn't malloc %d bytes", MAXNAMELEN) ;
     exit (-1) ;
     }
     
     if (!(seq = malloc (MALLOCBLOCK+1)))
       { fprintf (stderr, "Couldn't malloc %d bytes", MALLOCBLOCK) ;
	 exit (-1) ;
       }

  if (!(fil = fopen ( argv[1], "r" ))) 
    usage ();
  
  while ( readSequence(fil, dna2textConv, &seq, &seqname, &desc, &length) ) 
    /* once through per sequence */
    { 
      i = 0 ;
      while ( seqname[i] != ' ' && seqname[i] != '\0' && i < 256 )
	++i ;
      seqname[i] = '\0' ;

	  FILE *fp;

	  char *pathfile = "D:\\results\\";

	  char *outputFile = join1(pathfile, argv[1]);
		
	  if((fp=fopen(outputFile,"w"))==NULL){
		  printf("File cannot be opened\n");
	      exit (0);
	  }else{
	      printf("File opened for writing\n");
	  }		  
      fclose(fp);
	  if((fp=fopen(outputFile,"w+"))==NULL){
		  printf("File cannot be opened\n");
	      exit (0);
	  }else{
	      printf("File opened for writing\n");
	  }		  

      findspans (fp, 0, length, seq, seqname ) ;
	  fclose(fp);
    }

  exit (0);
}

void findspans (FILE *fp, int start, int end, char *seq, char *seqname )
{
  int i ;
  int sc = 0 ;
  int lsc = 0 ; 
  int imn = -1 ;  /* else sequences which start with pattern reported badly */
  int imx = 0 ;
  int mx = 0 ;
  int winlen = 0;
  float expect, obsToExp;

  int ncpg, ngpc, ngc, ng, nc;  
  int nacga, nacgc, nacgg,nacgt, nccga, nccgc, nccgg, nccgt, ngcga,ngcgc,ngcgg,ngcgt, ntcga,ntcgc,ntcgg, ntcgt;

  i = start ;
  while ( i < end )  
    {
      lsc = sc ;
      sc += ( end-1-i && seq[i]=='C' && seq[i+1]=='G' ) ? CPGSCORE : -1 ;
      sc = sc < 0 ? 0 : sc ;
/*      printf("%d \t %d \t%d \t %d \t%d \t%d\n", i, sc, lsc, imn, imx, mx) ; */
      if ( sc == 0 && lsc )  /* should threshold here */
	{
	  /* imn+1 is the start of the match. 
	     imx+1 is the end of the match, given pattern-length=2.
	     fetchdb using reported co-ordinates will return a sequence
	     which both starts and ends with a complete pattern.
	     Further +1 adjusts so that start of sequence = 1 
	  */

	  getstats (imn+1, imx+2, seq, seqname, &ncpg, &ngpc, &ng, &nc,&nacga, &nacgc, &nacgg,&nacgt, &nccga, &nccgc, &nccgg, &nccgt, &ngcga,&ngcgc,&ngcgg,&ngcgt, &ntcga,&ntcgc,&ntcgg, &ntcgt) ;
	  ngc = ng + nc;
      if (((imx+2)-(imn+2))>199 && (ngc*100.0/(imx+1-imn))>50.0 ) {
	/* old gos estimate	  printf("%s\t %d\t %d\t %d\t CpG: %d\t %.1f\t %.1f\n", seqname, imn+2, imx+2, mx, ncpg, ngc*100.0/(imx+1-imn), 1.0*ncpg/ngpc) ; */
	winlen=imx+1-imn;
	/* ASH 3/23/04: expected val from Gardiner-Garden & Frommer '87: */
	expect = (float)(nc * ng) / (float)winlen;
	obsToExp = (float)ncpg / expect;
			if ( obsToExp > 0.60 )
			{
			   printf("%s\t %d\t %d\t %d\t CpG: %d\t %.1f\t %.2f\t %.2f\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t\n",seqname, imn+2, imx+2, mx, ncpg, ngc*100.0/(imx+1-imn),1.0*ncpg/ngpc, obsToExp, nacga, nacgc, nacgg, nacgt, nccga, nccgc, nccgg, nccgt, ngcga, ngcgc, ngcgg, ngcgt, ntcga, ntcgc, ntcgg, ntcgt ) ; 
			   fprintf(fp,"%s\t %d\t %d\t %d\t CpG: %d\t %.1f\t %.2f\t %.2f\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t\n", seqname, imn+2, imx+2, mx, ncpg, ngc*100.0/(imx+1-imn), 1.0*ncpg/ngpc, obsToExp, nacga, nacgc, nacgg, nacgt, nccga, nccgc, nccgg, nccgt, ngcga, ngcgc, ngcgg, ngcgt, ntcga, ntcgc, ntcgg, ntcgt );
			}

      }
/* 	  printf("%s \t %d\t %d\t %d \n", seqname, imn+2, imx+2, mx ) ; 
  */
	  /* Recursive call searches from one position after the end of the 
	     last match to the current position */
	  findspans(fp,imx+2, i, seq, seqname ) ;
	  sc = lsc = imn = imx =  mx = 0 ;
	}
      imx = sc > mx ? i : imx ;
      mx = sc > mx ? sc : mx ;
      imn = sc == 0 ? i : imn ;
      ++i ;
    }
  if ( sc != 0 )  /* should threshold here */
    {
/*      printf("%d \t %d \t%d \t %d \t%d \t%d\n", i, sc, lsc, imn, imx, mx) ;  */
      
      /* ASH 3/23/04: Make this test & output consistent w/above. */
	  getstats (imn+1, imx+2, seq, seqname, &ncpg, &ngpc, &ng, &nc,&nacga, &nacgc, &nacgg,&nacgt, &nccga, &nccgc, &nccgg, &nccgt, &ngcga,&ngcgc,&ngcgg,&ngcgt, &ntcga,&ntcgc,&ntcgg, &ntcgt) ;
      ngc = nc + ng;
      if (((imx+2)-(imn+2))>199 && (ngc*100.0/(imx+1-imn))>50.0 ) {
	winlen=imx+1-imn;
	expect = (float)(nc * ng) / (float)winlen;
	obsToExp = (float)ncpg / expect;
			if ( obsToExp > 0.60 )
			{
				   printf("%s\t %d\t %d\t %d\t CpG: %d\t %.1f\t %.2f\t %.2f\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t\n", seqname, imn+2, imx+2, mx, ncpg, ngc*100.0/(imx+1-imn),1.0*ncpg/ngpc, obsToExp, nacga, nacgc, nacgg, nacgt, nccga, nccgc, nccgg, nccgt, ngcga, ngcgc, ngcgg, ngcgt, ntcga, ntcgc, ntcgg, ntcgt ) ; 
				   fprintf(fp,"%s\t %d\t %d\t %d\t CpG: %d\t %.1f\t %.2f\t %.2f\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t\n", seqname, imn+2, imx+2, mx, ncpg, ngc*100.0/(imx+1-imn), 1.0*ncpg/ngpc, obsToExp, nacga, nacgc, nacgg, nacgt, nccga, nccgc, nccgg, nccgt, ngcga, ngcgc, ngcgg, ngcgt, ntcga, ntcgc, ntcgg, ntcgt );
			}
      }
      
/*      printf("%s \t %d\t %d\t %d \n", seqname, imn+2, imx+2, mx ) ; 
   */
      findspans(fp,imx+2, end, seq, seqname ) ;
    }
}

/* ASH 3/23/04: Return separate counts of G's and C's for expected val calc. */
void getstats (int start, int end, char *seq, char *seqname, int *ncpg, int *ngpc, int *ng, int *nc, int *nacga, int *nacgc, int *nacgg, int *nacgt, int *nccga, int *nccgc, int *nccgg, int *nccgt, int *ngcga,int *ngcgc,int *ngcgg,int *ngcgt, int *ntcga,int *ntcgc,int *ntcgg, int *ntcgt )
{

int i ;
*ncpg = *ngpc = *ng = *nc = 0 ;
*nacga= *nacgc=*nacgg= *nacgt= *nccga= *nccgc= *nccgg= *nccgt= *ngcga= *ngcgc= *ngcgg= *ngcgt= *ntcga= *ntcgc= *ntcgg= *ntcgt =0;

for ( i = start ; i < end ; ++i ) 
  {

   if ( end-1-i && (seq[i]=='C' || seq[i]=='c') && (seq[i+1]=='G'||seq[i+1]=='g') ) ++*ncpg ; 
   if ( end-1-i &&  (seq[i]=='G'||seq[i]=='g')&& (seq[i+1]=='C' || seq[i+1]=='c') ) ++*ngpc ; 
   if ( seq[i]=='G' ) ++*ng ; 
   if ( seq[i]=='C' ) ++*nc ; 

      if(i>0 && i<end-1 && (seq[i]=='C'|| seq[i]=='c') && (seq[i+1]=='G'||seq[i+1]=='g') ){
		 char seqTemp[] = {"NNNN"};
		 seqTemp[0] = seq[i-1];
		 seqTemp[1] = seq[i];
		 seqTemp[2] = seq[i+1];
		 seqTemp[3] = seq[i+2];
        // char *seqTemp = strcat(seq[i-1],seq[i],seq[i+1],seq[i+2]);
	     if (strcasecmp(seqTemp, "ACGA") == 0) ++*nacga; 
		 if (strcasecmp(seqTemp, "ACGC") == 0) ++*nacgc; 
		 if (strcasecmp(seqTemp, "ACGG") == 0) ++*nacgg; 
		 if (strcasecmp(seqTemp, "ACGT") == 0) ++*nacgt; 
		 if (strcasecmp(seqTemp, "CCGA") == 0) ++*nccga; 
		 if (strcasecmp(seqTemp, "CCGC") == 0) ++*nccgc; 
		 if (strcasecmp(seqTemp, "CCGG") == 0) ++*nccgg; 
		 if (strcasecmp(seqTemp, "CCGT") == 0) ++*nccgt; 
		 if (strcasecmp(seqTemp, "GCGA") == 0) ++*ngcga; 
		 if (strcasecmp(seqTemp, "GCGC") == 0) ++*ngcgc; 
		 if (strcasecmp(seqTemp, "GCGG") == 0) ++*ngcgg; 
		 if (strcasecmp(seqTemp, "GCGT") == 0) ++*ngcgt; 
		 if (strcasecmp(seqTemp, "TCGA") == 0) ++*ntcga; 
		 if (strcasecmp(seqTemp, "TCGC") == 0) ++*ntcgc; 
		 if (strcasecmp(seqTemp, "TCGG") == 0) ++*ntcgg; 
		 if (strcasecmp(seqTemp, "TCGT") == 0) ++*ntcgt; 
	  }
   }
}

char *join1(char *a, char *b){
     char *c = (char *) malloc(strlen(a) + strlen(b) + 1); 
     if(c==NULL) exit(1);
     char *tempc =c;
     while(*a!='\0'){
         *c++ = *a++;
      }
     while((*c++=*b++)!='\0'){
           ;
      }
      return tempc;//
}