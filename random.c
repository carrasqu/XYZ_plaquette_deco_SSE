/***************************************************************************

	    LAGGED FIBONNACCI PSEUDO-RANDOM NUMBER GENERATOR
	    ================================================


                            Paul Coddington,
                 Northeast Parallel Architectures Center,
              Syracuse University, Syracuse NY 13244, U.S.A.
                           paulc@npac.syr.edu

                              October 1993.


BACKGROUND
----------

This program implements a lagged Fibonnacci pseudo-random number generator
using multiplication, with lags p and q, i.e. F(p,q,*) in the notation of
[Marsaglia].

A sequence of odd random integers X  is obtained from the formula
                                   n

   X  = X    * X        mod M
    n    n-p    n-q

where M is a large integer less than the machine precision (M is taken
to be 2^{31} here).

Unlike lagged Fibonnacci generators which use subtraction, addition, or
XOR, this generator seems to give "good" random numbers for ANY lag, not
necessarily large (see [Coddington]). However larger lags give larger
periods and empirically seem to give "better" random numbers (see
[Marsaglia], [Knuth], [Coddington]).

The period is given by (2^{p}-1)*2^{b-3}, where p is the largest lag
and b is the number of bits in the integers X_{n}, i.e. b = log_{2}(M).
The lag needs to be chosen large enough to give a period much longer
than the number of random numbers to be generated.

Only certain lags (p,q) will give the maximal period of the generator
(see [Marsaglia], [Knuth]). Here are some of those lags, with the
corresponding periods for 32-bit integers. As a comparison, note that
linear congruential generators have a period of about 10^{9} (which is
too short for most purposes), or 10^{18} for 64-bit integers (which is
fine for most purposes). A Teraflop-year is about 10^{21} floating point
operations.

    P     Q        period
  1279   418       10^{393}
   607   273       10^{191}
   521   168       10^{166}
   250   103       10^{84}
   127    63       10^{46}
    97    33       10^{38}
    89    38       10^{35}
    55    24       10^{25}
    43    22       10^{21}
    31    13       10^{18}
    24    10       10^{16}
    17     5       10^{14}
     7     3       10^{11}
     5     2       10^{10}

I would recommend using at least (43,22), and ideally (607,273). This
will only add an extra 607x4 = 2428 bytes to the memory required for
the user program, which should be completely negligible.

Note that very little is known theoretically about this generator,
however it performs very well in all empirical statistical tests (see
[Marsaglia], [Coddington]).
THIS DOES NOT NECESSARILY MEAN IT WILL PERFORM WELL FOR YOUR APPLICATION.
It is a good idea to check by using another good generator such as RANF
or DRAND48 or a good 64-bit generator (see [Coddington]).



THE PROGRAM
-----------

This program is based on the implementation of RANMAR in [James].
It is the same as RANMAR except it uses multiplication instead of
subtraction, odd integers instead of reals, and does not include an
extra Weyl generator.

The program accesses the "lag table", i.e. the array storing the previous
p values in the sequence, using a "circular list", i.e. the elements of
the array are not moved with every new addition, only the pointers pt0
and pt1 to the (n-p)th and (n-q)th value are changed.

NOTE THAT YOU *MUST* CALL THE INITIALIZATION SUBROUTINE (rand_init) BEFORE
USING THE GENERATOR.

In order for the sequence to be random, the initial values of the lag
table have to be "random enough". There have not to my knowledge been any
quantitive tests of this point, however initializing the lag table using
a good linear congruential generator (the Unix/C generator RAND is used
here) seems empirically to work fine.

The computation required to generate a new random number, apart from
updating 2 pointers, is 1 integer multiplication, which multiplies two
32-bit integers and returns the answer modulo 2^{32} (this is the standard
for unsigned ints in C), and 1 floating point multiplication, to convert
the integers into a floating point number in [0,1).

For simplicity and clarity the lag table is indexed from 1 to P rather
than 0 to P-1.



REFERENCES
----------

For more details, see:


P.D. Coddington, "Monte Carlo tests of random number generators",
NPAC technical report SCCS-526.

F. James, "A Review of Pseudo-random Number Generators",
Comput. Phys. Comm. 60 (1990) 329.

D.E. Knuth, "The Art of Computer Programming Vol. 2: Seminumerical
Methods", (Addison-Wesley, Reading, Mass., 1981).

P. L'Ecuyer, "Random numbers for simulation", Comm. ACM 33:10 (1990) 85.

G.A. Marsaglia, "A current view of random number generators",
in Computational Science and Statistics: The Interface,
ed. L. Balliard (Elsevier, Amsterdam, 1985).


***************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#define TWOTO31 2147483648.0
double  TWOTONEG31=(1.0/TWOTO31);

/* Lags */
#define P 1279
#define Q 418

/* Variables for lagged Fibonacci generator */
unsigned int u[P+1];
int pt0, pt1;

/* Variables for the linear congruential generator used in initialization */
#define M ( 0x7fffffff )     /* 2^{31}-1 -- mask by M to get modulo 2^{31} */
#define	A 1103515245
#define	C 12345


/***************************************************************************
   Initialize the random number generator.

   Every bit of all the seeds in the lag table is initialized using the
   linear congruential generator RAND (the standard Unix/C 32-bit generator)

        s    = (A * s  + C)   mod M
         i+1         i

   A = 1103515245   C = 12345   M = 2^{31}-1

   This routine MUST be called once at the beginning of the program
   before the random number generator is called.
***************************************************************************/
#ifdef RISC
  void rand_init(seed)
#else
  void rand_init_(seed)
#endif
int * seed;
{
   int i,j;
   int s,t,randx;

   /* Initialize the LCG RAND */
   randx = *seed;
   /* printf("%ld\n",*seed) ; */

   /* Initialize every bit of the lag table using RAND */
   for (i=1;i<=P;i++) {
      s = 0;
      t = 1;
      for (j=0;j<32;j++) {
         randx = (A * randx + C) & M;  /* RAND */
         if ( ( (double)randx * TWOTONEG31 ) < 0.5 ) {
            s = s + t;
         }
	 t = t << 1;
      }
      /* Seeds must be odd integers, so make the last bit a 1 */
      s = s | 1;
      u[i] = s;
   }

   /* Initialize pointers */
   pt0 = P;
   pt1 = Q;

   return;
}

/***************************************************************************
   Return a random double in [0.0, 1.0).

   Note that only the most significant 31 bits of drand1 will actually
   be random, since a random 32-bit integer is converted into a double.

   To return a float instead of a double, simply change the type of rand
   and drand1 to float.
***************************************************************************/
#ifdef RISC
  double drand1()
#else
  double drand1_()
#endif
{
   unsigned int uni;
   double rand;

   /*
    * Note that u[i] and uni are unsigned ints so multiplication
    * modulo 2^{32} is done automatically in C.
    */
   uni = u[pt0] * u[pt1];
   u[pt0] = uni;
   pt0 --;
   if (pt0 == 0) pt0 = P;
   pt1 --;
   if (pt1 == 0) pt1 = P;
   /*
    * Shift to get rid of last odd bit.
    */
   uni = uni >> 1;
   /*
    * uni is now between 0 and 2^{31-1}, so multiply by 2^{-31}
    * to give a random floating point number in [0,1).
    */
   rand = (double)uni * TWOTONEG31;
   return(rand);
}

/***************************************************************************
   Store the lag table and pointers.
***************************************************************************/
write_seeds_(seedfilename)
char seedfilename[128];
{
   FILE *f1,*fopen();
   int i,pointer[2];

   pointer[0] = pt0;
   pointer[1] = pt1;

   /* Open a file to store the lag table and pointers */
   f1 = fopen(seedfilename,"w") ;
   for (i=1;i<=P;i++) {
      fprintf(f1,"%lu\n",u[i]);
   }
   fflush(f1);
   fprintf(f1,"%d %d\n",pointer[0],pointer[1]);
   fclose(f1);

   return 0;
}

/***************************************************************************
   Read in stored lag table and pointers.
***************************************************************************/
read_seeds_(seedfilename)
char seedfilename[128];
{
   FILE *f1,*fopen();
   int i,pointer[2];

   /* Open the file containing the lag table and pointers */
   f1 = fopen(seedfilename,"r");
   for (i=1;i<=P;i++) {
      fscanf(f1,"%lu",&u[i]);
   }
   fscanf(f1,"%d %d",&pointer[0],&pointer[1]);
   fclose(f1);

   pt0 = pointer[0];
   pt1 = pointer[1];

   return 0;
}


