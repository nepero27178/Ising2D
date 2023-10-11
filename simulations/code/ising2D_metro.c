#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<time.h>

#include"./random.h"

#define SIZE 4              // Square lattice SIZE x SIZE
#define STRING_LENGTH 50    // Maximum lenght of datafile name

// ----- Functions -----

/* 
@ Total magnetization calculator function:
This function spans all over the lattice and sums the value of each site to the
total magnetization "sum".
*/

double magn(int lattice[SIZE][SIZE])
  {
  long int rx, ry, sum;

  sum=0;
  for(rx=0; rx<SIZE; rx++)
     {
     for(ry=0; ry<SIZE; ry++)
        {
        sum+=lattice[rx][ry];
        }
     }

  return (double) sum / (double) (SIZE*SIZE);
  }


/* 
@ Total energy calculator function:
This function spans all over the lattice and sums the value of each site to the
total magnetization "sum". Note that the sum is only performed on (rx,ry),
(rx+1,ry) and (rx,ry+1). This is because che contribution of (rx-1,ry) and
(rx,ry-1) was already counted in the previous steps.
*/

double energy(int lattice[SIZE][SIZE]) 
  {
  long int rx, ry, tmp, sum;

  sum=0;
  for(rx=0; rx<SIZE; rx++)
     {
     for(ry=0; ry<SIZE; ry++)
        {
        tmp=(rx+1)%SIZE;
        sum+=-lattice[rx][ry]*lattice[tmp][ry];
       
        tmp=(ry+1)%SIZE;
        sum+=-lattice[rx][ry]*lattice[rx][tmp];
        }
     }

  return (double) sum / (double) (SIZE*SIZE);
  }


/*
@ Metropolis steps:
1) Evaluate the exchange energy with the four neighbours (store in "sumnn");
2) If the central spin is flipped with respect to the total neighbouring spin,
then the flip is physically favoured because it reduces the lattice energy,
therefore the step is automatically accepted;
3) In the other case, we must sort a random number to accept the step with
respect to the ratio of boltzmann probabilities of the configurations. 
*/

int metropolis(int lattice[SIZE][SIZE], 
               long int rx,
               long int ry,
               double const * const restrict acc_prob)
  {
  long int tmp;
  int sumnn, acc=0;

  // Start with sumnn=0
  sumnn=0;
  
  // Sum all the neighbours
  tmp=(rx+1)%SIZE;
  sumnn+=lattice[tmp][ry];

  tmp=(rx-1);
  if(tmp<0)
    {
    tmp=SIZE-1;
    }
  sumnn+=lattice[tmp][ry];

  tmp=(ry+1)%SIZE;
  sumnn+=lattice[rx][tmp];

  tmp=(ry-1);
  if(tmp<0)
    {
    tmp=SIZE-1;
    }
  sumnn+=lattice[rx][tmp];
  
  // Multiply by the center
  sumnn*=lattice[rx][ry];

  // Metropolis step
  if(sumnn<0)
    {
    lattice[rx][ry]=-lattice[rx][ry];
    acc=1;
    }
  else
    {
    // acc_prob matrix pre-calclulated (see "main")
    if(myrand()<acc_prob[sumnn])
      {
      lattice[rx][ry]=-lattice[rx][ry];
      acc=1; // Update acceptance
      }
    }
  
  // Exit: update lattice and acceptance counter
  return acc;
  }


// ----- Main -----

int main(int argc, char **argv)
    {
    int i, lattice[SIZE][SIZE];
    long int rx, ry, rxaux, ryaux, sample, iter, acc; 
    double beta, locE, locM;
    double acc_prob[5];
  
    char datafile[STRING_LENGTH];
    FILE *fp;
    
    // Seeds to pass to random generator: use time
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 4)
      // Help message
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s beta sample datafile\n\n", argv[0]);
      fprintf(stdout, "  beta = inverse temperature\n");
      fprintf(stdout, "  sample = number of drawn to be extracted\n");
      fprintf(stdout, "  datafile = name of the file on which to write the data\n");
      fprintf(stdout, "  size of the lattice define by macro: now it is %d\n\n", SIZE);
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  E, M (E=energy per site, M=magnetization per site), one line for each draw\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // Read input values 
      beta=atof(argv[1]);
      sample=atol(argv[2]);
        
      // Read datafile name
      if(strlen(argv[3]) >= STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      else
        {
        // Copy the input to the variable "datafile"
        strcpy(datafile, argv[3]);
        }
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // Initialize random number generator
    myrand_init(seed1, seed2);

    // Initialize lattice to ordered start
    for(rx=0; rx<SIZE; rx++)
       {
       for(ry=0; ry<SIZE; ry++)
          {
          // Start with all spin "up"
          lattice[rx][ry]=1;
          }
       }

    // Open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    /* 
    @ Initialize acceptance probability:
    We have 10 possible configuration. Using Z2 simmetry we reduce those to 9.
    All the configuration with negative difference in exchange energy are
    automatically accepted, so we are left with those where the flip goes
    agains the mean trend of the neighbours.
    Since the difference in energy can vary discretely, we precalculate all
    the possibilities and store them in "acc_prob".
    */    
    for(i=0; i<5; i++)
       {
       acc_prob[i]=exp(-2.0*beta*((double)i));
       }
    
    /*
    @ Main run:
    1) Repeat sampling "sample" number of times;
    2) For each sampling span all over the lattice and for each site extract
    another site at random on the lattice (rxaux,ryaux);
    3) Run here the metropolis algorithm;
    4) Once the spanning over sites is complete, L^2 steps have been done.
    Calculate energy and magnetization and store them.
    */
    acc=0;
    for(iter=0; iter<sample; iter++)
       {
       for(rx=0; rx<SIZE; rx++)
          {
          for(ry=0; ry<SIZE; ry++)
             {
             rxaux=(long int)((double)SIZE * myrand());
             ryaux=(long int)((double)SIZE * myrand());
             acc+=metropolis(lattice, rxaux, ryaux, acc_prob);
             }
          }

       locE=energy(lattice);
       locM=magn(lattice);

       fprintf(fp, "%.12f %.12f\n", locE, locM);
       }

    // Close datafile
    fclose(fp);


    /*
    Acceptance rate: mean number of acceptance probability per site
    Therefore: AR = "number of accepted steps" over "total number of spins 
    spanned" which is "number of samples" * "number of spins per lattice")    
    */
    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) (SIZE*SIZE));

    return EXIT_SUCCESS;
    }


