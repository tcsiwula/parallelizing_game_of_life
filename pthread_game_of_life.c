 /*=============================================================================
 File:         =        pthread_game_of_life.c
 Created:      =        12/08/2015
 Author:       =        Tim Siwula <tcsiwula@usfca.edu>
 University:   =        University of San Francisco
 Class:        =        Computer Science 220: Introduction to Parallel Computing
 Liscense:     =        GPLv2
 Version:      =        0.03
 ==============================================================================
 Purpose:      =        Implement John Conway's Game of Life.  The game is
    `    `played'' on a board or *world* consisting of a rectangular grid with
 *           m rows and n columns.  Each cell in the grid is ``alive''
 *           or ``dead.'' An initial generation (generation 0) is either
 *           entered by the user or generated using a random number
 *           generator.
 *
 *           Subsequent generations are computed according to the
 *           following rules:
 *
 *              - Any live cell with fewer than two live neighbors
 *              dies, as if caused by under-population.
 *              - Any live cell with two or three live neighbors
 *              lives on to the next generation.
 *              - Any live cell with more than three live neighbors
 *              dies, as if by over-population.
 *              - Any dead cell with exactly three live neighbors
 *              becomes a live cell, as if by reproduction.
 *
 *           Updates take place all at once.
 ==============================================================================
 Compile:      =        gcc -g -Wall -o life pthread_game_of_life.c
 ==============================================================================
 Debug:        =        gcc -DDEBUG -g -Wall -o life pthread_game_of_life.c
 ==============================================================================
 Run:          =        ./life <r> <s> <m> <n> <max> <'i'|'g'> -lpthread
 ==============================================================================
 Run args:     =        1. r: number of rows of threads
               =        2. s: number of columns of threads
               =        3. m: number of rows in the world
               =        4. n: number columns in world
               =        5. max: limit for generations
               =        6. <'i'|'g'>: g = generate, r = read initial generation.
 ==============================================================================
 Input:        =        If command line has the "input" char ('i'), the first
 *                      generation.  Each row should be entered on a separate
 *                      line of input.  Live cells should be indicated with
 *                      a capital 'X', and dead cells with a blank, ' '.
 *                      If command line had the "generate" char ('g'), the
 *                      probability that a cell will be alive.
 ==============================================================================
 Output:       =        The initial world (generation 0) and the world after
               =        each subsequent generation up to and including
                        generation = max_gen.  If all of the cells die,
                        the program will terminate.
 ==============================================================================
 Note:         =        1. We use a one dimensional array to store the worlds.
               =        2. This implementation uses a "toroidal world" in which
                           the last row of cells is adjacent to the first row,
                           and the last column of cells is adjacent to the first
 ............................................................................ */

/* ========================================================================== */
/*                          External libaries                                 */
/* ========================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
/* .......................................................................... */

/* ========================================================================== */
/*                          Constants                                         */
/* ========================================================================== */
#define LIVE 1
#define DEAD 0
#define LIVE_IO 'X'
#define DEAD_IO ' '
#define MAX_TITLE 1000
/* .......................................................................... */

/* ========================================================================== */
/*                          Custom function definitions                       */
/* ========================================================================== */
void Usage(char prog_name[]);
void Read_world(char prompt[], int wp[], int m, int n);
void Gen_world(char prompt[], int wp[], int m, int n);
void Print_world(char title[], int wp[], int m, int n);
void* Play_life(void* rank);
void Synchronize(int live_count, int curr_gen);
int  Count_nbhrs(int *wp, int m, int n, int i, int j);
void Get_args(int argc, char* argv[]);
/* .......................................................................... */

/* ========================================================================== */
/*                          Shared-Global variables                       */
/* ========================================================================== */
int r, s, m, n, max_gens;
char ig;
int thread_count, rowSize, colSize;
int *w1, *w2, *wp, *twp;
int global_live_count = 0, barrier_thread_count = 0;//, curr_gen = 0;
pthread_mutex_t barrier_mutex;
pthread_cond_t ok_to_proceed;
/* .......................................................................... */

/* ========================================================================== */
/*                          Main( )                                           */
/* ========================================================================== */
int main(int argc, char* argv[]) {
   Get_args(argc, argv);
   long thread;
   pthread_t* thread_handles = malloc(thread_count*sizeof(pthread_t));
   pthread_mutex_init(&barrier_mutex, NULL);
   pthread_cond_init(&ok_to_proceed, NULL);
   w1 = malloc(m*n*sizeof(int));
   w2 = malloc(m*n*sizeof(int));
   wp = w1;
   twp = w2;

   if (ig == 'i')
      Read_world("Enter generation 0\n", wp, m, n);
   else
      Gen_world("What's the prob that a cell is alive?\n", wp, m, n);

    for (thread = 0; thread < thread_count; thread++)
          pthread_create(&thread_handles[thread], (pthread_attr_t*) NULL,
                            Play_life, (void*) thread);

    for (thread = 0; thread < thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);

   pthread_mutex_destroy(&barrier_mutex);
   pthread_cond_destroy(&ok_to_proceed);
   free(thread_handles);
   free(w1);
   free(w2);
   return 0;
}  /* main */
/* .......................................................................... */

/* ========================================================================== */
/*                          Play_life( )                                      */
/* ========================================================================== */
 /* Function:     Play_life
 * Purpose:      Play Conway's game of life.  (See header doc)
 * In args:      rank:  a threads rank
 */
void* Play_life(void* rank)
{
   int i, j, count, live_count, curr_gen = 0;
   int startRow, startCol, stopRow, stopCol;
   long my_rank = (long) rank;
   char title[MAX_TITLE];

   if(my_rank == 0)
   {
     sprintf(title, "Generation %d of %d", curr_gen, max_gens);
     Print_world(title, wp, m, n);
   }
       curr_gen++;


   rowSize = m/r;
   colSize = n/s;

   startRow = (my_rank % s) * rowSize;
   startCol = (my_rank % s) * colSize;

   stopRow = rowSize + r;
   stopCol = colSize + s;


   while (curr_gen <= max_gens)
   {
      live_count = 0;
      for (i = startRow; i < stopRow; i++)
      {
         for (j = startCol; j < stopCol; j++)
         {
            count = Count_nbhrs(wp, m, n, i, j);
        if (count < 2 || count > 3)
           twp[i*n + j] = DEAD;
        else if (count == 2)
           twp[i*n + j] = wp[i*n + j];
        else /* count == 3 */
           twp[i*n + j] = LIVE;
        if (twp[i*n + j] == LIVE) live_count++;
     }
  }
        Synchronize(live_count, curr_gen);
      curr_gen++;
     if (global_live_count == 0) break;

   }
   return 0;
}  /* Play_life */
/* .......................................................................... */

/* ============================================================================
Function:     =        Synchronize()
==============================================================================
Purpose:      =        Uses conditional barries to sync the threads.
==============================================================================
Input arg:    =        1.live_count: A given threads current gen alive count.
              =        2. curr_gen: Which generation the threads are on.
=========================================================================== */
void Synchronize(int live_count, int curr_gen)
{
   char title[MAX_TITLE];
   int *tmp;
   pthread_mutex_lock(&barrier_mutex);
   barrier_thread_count++;
   global_live_count += live_count;

  if (barrier_thread_count == thread_count)
  {
     tmp = wp;
     wp = twp;
     twp = tmp;

     if (global_live_count != 0 && live_count != 0) {
      sprintf(title, "Generation %d of %d", curr_gen, max_gens);
      Print_world(title, wp, m, n);
     }
          barrier_thread_count = 0;

     pthread_cond_broadcast(&ok_to_proceed);
  } else
  {
     // Wait unlocks mutex and puts thread to sleep.
     //    Put wait in while loop in case some other
     // event awakens thread.
     while (pthread_cond_wait(&ok_to_proceed,
               &barrier_mutex) != 0);
     // Mutex is relocked at this point.
  }
  pthread_mutex_unlock(&barrier_mutex);
}
/* .......................................................................... */

/* ============================================================================
 Function:     =        Get_args
 ==============================================================================
 Purpose:      =        Get the command line arguments.
 ==============================================================================
 Input arg:  =        1. argc: Flag for list creation type.
             =        2. *argv[]: A pointer to an array.
             =        5. list_size: The size of the processors array.
             =        6. number_of_processors: The size of the processors array.
 ==============================================================================
 Out args:   =        1. list_size_p: A pointer an integer.
             =        2. generate_list_p: pointer to an i or g for list type.
 =========================================================================== */
void Get_args(int argc, char* argv[])
{
     if (argc > 8 || argc < 7) Usage(argv[0]);
     r = strtol(argv[1], NULL, 10);
     s = strtol(argv[2], NULL, 10);
     m = strtol(argv[3], NULL, 10);
     n = strtol(argv[4], NULL, 10);
     max_gens = strtol(argv[5], NULL, 10);
     ig = argv[6][0];
     thread_count = r*s;
}  /* Get_args */
/* .......................................................................... */

/*---------------------------------------------------------------------
 * Function:   Usage
 * Purpose:    Show user how to start the program and quit
 * In arg:     prog_name
 */
void Usage(char prog_name[]) {
   fprintf(stderr, "usage: %s <rows> <cols> <max> <i|g>\n", prog_name);
   fprintf(stderr, "    r = number of rows of threads\n");
   fprintf(stderr, "    s = number of columns of threads\n");
   fprintf(stderr, "    m = number of rows in the world\n");
   fprintf(stderr, "    n = number of cols in the world\n");
   fprintf(stderr, "     max = max number of generations\n");
   fprintf(stderr, "       i = user will enter generation 0\n");
   fprintf(stderr, "       g = program should generate generation 0\n");
   exit(0);
}  /* Usage */

/*---------------------------------------------------------------------
 * Function:   Read_world
 * Purpose:    Get generation 0 from the user
 * In args:    prompt
 *             m:  number of rows in visible world
 *             n:  number of cols in visible world
 * Out arg:    wp:  stores generation 0
 *
 */
void Read_world(char prompt[], int wp[], int m, int n) {
   int i, j;
   char c;

   printf("%s\n", prompt);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         scanf("%c", &c);
         if (c == LIVE_IO)
            wp[i*n + j] = LIVE;
         else
            wp[i*n + j] = DEAD;
      }
      /* Read end of line character */
      scanf("%c", &c);
  }
}  /* Read_world */


/*---------------------------------------------------------------------
 * Function:   Gen_world
 * Purpose:    Use a random number generator to create generation 0
 * In args:    prompt
 *             m:  number of rows in visible world
 *             n:  number of cols in visible world
 * Out arg:    wp:  stores generation 0
 *
 */
void Gen_world(char prompt[], int wp[], int m, int n) {
   int i, j;
   double prob;
#  ifdef DEBUG
   int live_count = 0;
#  endif

   printf("%s\n", prompt);
   scanf("%lf", &prob);

   srandom(1);
   for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
         if (random()/((double) RAND_MAX) <= prob) {
            wp[i*n + j] = LIVE;
#           ifdef DEBUG
            live_count++;
#           endif
         } else {
            wp[i*n + j] = DEAD;
         }

#  ifdef DEBUG
   printf("Live count = %d, request prob = %f, actual prob = %f\n",
         live_count, prob, ((double) live_count)/(m*n));
#  endif
}  /* Gen_world */


/*---------------------------------------------------------------------
 * Function:   Print_world
 * Purpose:    Print the current world
 * In args:    title
 *             m:  number of rows in visible world
 *             n:  number of cols in visible world
 *             wp:  current gen
 *
 */
void Print_world(char title[], int wp[], int m, int n) {
   int i, j;
   printf("%s:", title);

   for (i = -1; i < m+1; i++)
   {
      for (j = -1; j < n+1; j++)
      {
             if(i >= 0 && i < m && j >= 0 && j < n)
             {
                 if (wp[i*n + j] == LIVE)
                    printf("%c", LIVE_IO);
                 else
                    printf("%c", DEAD_IO);
            }
      }
      printf("\n");
   }
}  /* Print_world */

/*---------------------------------------------------------------------
 * Function:   Count_nbhrs
 * Purpose:    Count the number of living nbhrs of the cell (i,j)
 * In args:    wp:  current world
 *             m:   number of rows in world
 *             n:   number of cols in world
 *             i:   row number of current cell
 *             j:   col number of current cell
 * Ret val:    The number of neighboring cells with living neighbors
 *
 * Note:       Since the top row of cells is adjacent to the bottom
 *             row, and since the left col of cells is adjacent to the
 *             right col, in a very small world, it's possible to
 *             count a cell as a neighbor twice.  So we assume that
 *             m and n are at least 3.
 */
int Count_nbhrs(int* wp, int m, int n, int i, int j) {
   int i1, j1, i2, j2;
   int count = 0;

   for (i1 = i-1; i1 <= i+1; i1++)
      for (j1 = j-1; j1 <= j+1; j1++) {
         i2 = (i1 + m) % m;
         j2 = (j1 + n) % n;
         count += wp[i2*n + j2];
      }
   count -= wp[i*n + j];

   return count;
}  /* Count_nbhrs */
