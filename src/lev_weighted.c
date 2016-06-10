#include<string.h>
#include<stdio.h>
#include <R.h>

/*This is a weighted levenshtein distance, whereby the 
  costs of symbol-symbol matches are provided
*/


/*
here, subst_cost is a (m+1)x(n+1) length vector specifying the
  relative distances between symbols.

the last row and column indicate the insertion costs--what happens
when a symbol is compared to 'NULL' essentially.

*/
int lookup(int row, int col, int maxi, int maxj);

double levenshtein_internal(int m, int n,
			    double* ins1,double * ins2,
			    double * sub_c, int print);

  
void lev_weighted(int *length_1, int *length_2,
		  double *sub_cost,
		  double *ins1, double* ins2,
		  int *print, double * ans);


void lev_weighted(int  *length_1, int  *length_2,
		  double * subst_cost,
		  double *ins1, double* ins2,
		  int *print,double * ans)
{


  //  Rprintf("lengths: %d, %d\n",*length_1, *length_2);
  //  Rprintf("inscost: %f\n",*ins_cost);

  int l = (*length_1) * (*length_2);

  int max_length= (*length_1) > (*length_2) ?(*length_1) : (*length_2);  
     
  double lev_dist=levenshtein_internal(*length_1,*length_2,
				       ins1,ins2, subst_cost, *print);
  *ans =lev_dist;
  // 		Rprintf("Vergleiche %s, %s\n",str_1, str_2); // Debug-Ausgabe
  // 		Rprintf("Levenshtein-Distanz: %d\n",lev_dist);   
  // 		Rprintf("Levenshtein-Distanz: %d\n",ans[str_ind]);   
} 

/*
 * Below follows extract of PostgreSQL modul fuzzystrmatch. Only the relevant
 * function levenshtein_internal is retained here. Some changes were made
 * to fit the needs of an R function (R_alloc is used, PostgreSQL 
 * error handling has been removed, static keyword is removed from signature).   
 */
 
/* Definition fom PostgreSQL sources */ 

#define Min(x, y)       ((x) < (y) ? (x) : (y))


/*
 * fuzzystrmatch.c
 *
 * Functions for "fuzzy" comparison of strings
 *
 * Joe Conway <mail@joeconway.com>
 *
 * $PostgreSQL: pgsql/contrib/fuzzystrmatch/fuzzystrmatch.c,v 1.32 2010/01/02 16:57:32 momjian Exp $
 * Copyright (c) 2001-2010, PostgreSQL Global Development Group
 * ALL RIGHTS RESERVED;
 *
 * levenshtein()
 * -------------
 * Written based on a description of the algorithm by Michael Gilleland
 * found at http://www.merriampark.com/ld.htm
 * Also looked at levenshtein.c in the PHP 4.0.6 distribution for
 * inspiration.
 * Configurable penalty costs extension is introduced by Volkan
 * YAZICI <volkan.yazici@gmail.com>.
 *
 * (comments on metaphone omitted)
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without a written agreement
 * is hereby granted, provided that the above copyright notice and this
 * paragraph and the following two paragraphs appear in all copies.
 *
 * IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
 * DOCUMENTATION, EVEN IF THE AUTHOR OR DISTRIBUTORS HAVE BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
 * ON AN "AS IS" BASIS, AND THE AUTHOR AND DISTRIBUTORS HAS NO OBLIGATIONS TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 *
 */

/*
 */
double
levenshtein_internal(int m, int n,
		     double* ins1,double * ins2,
                     double * sub_c,int print)

{
    double        *prev;
    double        *curr;
    int         i,
                j;

    //int print = 1;

	    
    const char *x;
    const char *y;

    if(print)
      {
	for(i=0; i<m; i++) Rprintf("%f ",ins1[i]);
	Rprintf("\n");

	for(i=0; i<n; i++) Rprintf("%f ",ins2[i]);
	Rprintf("\n");
      }
    
    /* One more cell for initialization column and row. */
    ++m;
    ++n;

    /*
     * Instead of building an (m+1)x(n+1) array, we'll use two different
     * arrays of size m+1 for storing accumulated values. At each step one
     * represents the "previous" row and one is the "current" row of the
     * notional large array.
     */
//    prev = (int *) palloc(2 * m * sizeof(int));
    prev = (double *) R_alloc(sizeof(double), 2 * (n+1));
    curr = prev + (n+1);

    if(print) Rprintf("----------------\n");
    prev[0] = 0;
    if(print) Rprintf("%1.3f ",prev[0]);

    /* Initialize the "previous" row to 0..cols */
    for (i = 1; i < n; i++)
      {
        prev[i] = prev[i-1] + ins2[i-1];
	if(print)Rprintf("%1.3f ",prev[i]);
      }

    if(print)Rprintf("\n");
    int counter = 0;

    /* Loop through rows of the notional array */
    for (j = 1; j < m; j++)
    {
        double *temp;


    /*
         * First cell must increment sequentially,
	 as we're on the j'th row of
         * the (m+1)x(n+1) array.
         */
        curr[0] = prev[0] +  ins1[j-1];
	if(print)Rprintf("%1.3f ",curr[0]);

        for (i = 1; i < n; i++)
        {
	  double         ins;  /*down*/
	  double         del;  /*right*/
	  double         sub;  /*diagonal*/

	  /* Calculate costs for probable operations. */
	  
	  /* Insertion; use insert/delete cost in matrix;
	     a move from above
	   */
	  
	  ins = prev[i] + ins1[j-1];

	  /* Deletion; use insertion/deletion cost in matrix;
	   a move from left*/
	  del = curr[i - 1] + ins2[i-1];

	  counter =lookup(j-1,i-1,m-1,n-1);
	  //Rprintf("count:  %d,%d=%d\n",j-1,i-1,counter);
	  
	  /*push; the current mismatch plus running mismatch*/
	  sub = prev[i - 1] + sub_c[lookup(j-1,i-1,m-1,n-1)];

	  
	  /*	  printf("internal5 %d %d %d:   %f %f %f-%f\n",
		  j,i,counter, ins,del,sub,sub_c[lookup(j-1,i-1,m-1,n-1)]);*/
	  
	    /*((*x == *y) ? 0 : sub_c);       /* Substitution */

            /* Take the one with minimum cost. */
            curr[i] = Min(ins, del);
            curr[i] = Min(curr[i], sub);
	    if(print) Rprintf("%1.3f ",curr[i]);
        }
	if(print)Rprintf("\n");
        /* Swap current row with previous row. */
        temp = curr;
        curr = prev;
        prev = temp;
    }

    /*
     * Because the final value was swapped from the previous row to the
     * current row, that's where we'll find it.
     */

    return prev[n - 1];
}


/*
  this looks up a value in the matrix, which has been
  flattened before sending to c compiled code.
  i and j should be 0-based.  maxi,maxj are 
  the size of the matrix (not 0-based) 10 is 10.
*/
int lookup(int row, int col, int maxi, int maxj)
{
  /*  Rprintf("lookup: %d %d %d %d \n",row,col,maxi,maxj);*/
  return  row * (maxj) + col;
  
}
