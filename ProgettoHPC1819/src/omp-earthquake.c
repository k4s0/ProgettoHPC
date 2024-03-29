/****************************************************************************
 *
 * earthquake.c - Simple 2D earthquake model
 *
 * Copyright (C) 2018 Moreno Marzolla <moreno.marzolla(at)unibo.it>
 * Last updated on 2018-12-29
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------------
 *
 * Versione di riferimento del progetto di High Performance Computing
 * 2018/2019, corso di laurea in Ingegneria e Scienze Informatiche,
 * Universita' di Bologna. Per una descrizione del modello si vedano
 * le specifiche sulla pagina del corso:
 *
 * http://moreno.marzolla.name/teaching/HPC/
 *
 * Per compilare:
 *
 * gcc -D_XOPEN_SOURCE=600 -std=c99 -Wall -Wpedantic earthquake.c -o earthquake
 *
 * (il flag -D_XOPEN_SOURCE=600 e' superfluo perche' viene settato
 * nell'header "hpc.h", ma definirlo tramite la riga di comando fa si'
 * che il programma compili correttamente anche se inavvertitamente
 * non si include "hpc.h", o per errore non lo si include come primo
 * file come necessario).
 *
 * Per eseguire il programma si puo' usare la riga di comando seguente:
 *
 * ./earthquake 100000 256 > out
 *
 * Il primo parametro indica il numero di timestep, e il secondo la
 * dimensione (lato) del dominio. L'output consiste in coppie di
 * valori numerici (100000 in questo caso) il cui significato e'
 * spiegato nella specifica del progetto.
 *
 ****************************************************************************/
 /*Lorenzo Casini 27/09/2019*/
#include "hpc.h"

#include <stdio.h>

#include <stdlib.h>     /* rand() */

#include <assert.h>

/* energia massima */
#define EMAX 4.0f
/* energia da aggiungere ad ogni timestep */
#define EDELTA 1e-4

#define HALO 1 /*ghost area*/
/**
 * Restituisce un puntatore all'elemento di coordinate (i,j) del
 * dominio grid con n colonne.
 */
static inline float *IDX(float * grid, int i, int j, int n) {
  return (grid + (i + HALO) * n + j + HALO);
}

/**
 * Restituisce un numero reale pseudocasuale con probabilita' uniforme
 * nell'intervallo [a, b], con a < b.
 */
float randab(float a, float b) {
  return a + (b - a) * (rand() / (float) RAND_MAX);
}

/**
 * Inizializza il dominio grid di dimensioni n*n con valori di energia
 * scelti con probabilità uniforme nell'intervallo [fmin, fmax], con
 * fmin < fmax.
 *
 * NON PARALLELIZZARE QUESTA FUNZIONE: rand() non e' thread-safe,
 * qundi non va usata in blocchi paralleli OpenMP; inoltre la funzione
 * non si "comporta bene" con MPI (i dettagli non sono importanti, ma
 * posso spiegarli a chi e' interessato). Di conseguenza, questa
 * funzione va eseguita dalla CPU, e solo dal master (se si usa MPI).
 */
void setup(float *grid, int n, float fmin, float fmax) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      * IDX(grid, i, j, n) = randab(fmin, fmax);
    }
  }
}

/**
 * Somma delta a tutte le celle del dominio grid di dimensioni
 * n*n. Questa funzione realizza il passo 1 descritto nella specifica
 * del progetto.
 * ------------------------------------------------------------------
 * Modifica apportata:
 * si presenta un caso di embarassingly parrallel structure, ho parallelizzato
 * utilizzando la clausola parrallel for di OpenMP e collassando 
 * i due cicli con aggiunta di clausola default(none) e specificando
 * lo scope delle variabili shared e private di default quelle dei due cicli. 
 */

void increment_energy(float *grid, int n, float delta) {
  #pragma omp parallel for collapse(2) default (none) shared(grid, n, delta)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      * IDX(grid, i, j, n) += delta;
    }
  }
}

/**
 * Restituisce il numero di celle la cui energia e' strettamente
 * maggiore di EMAX.
 * -------------------------------------------------------------
 * Modifica apportata:
 * parallelizzazione con collasso (2) e aggiunta di riduzione
 */
int count_cells(float *grid, int n) {
  int c = 0;
  #pragma omp parallel for collapse(2) default (none) shared(grid, n) reduction(+: c)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ( * IDX(grid, i, j, n) > EMAX) {
        c++;
      }
    }
  }
  return c;
}

/** 
 * Distribuisce l'energia di ogni cella a quelle adiacenti (se
 * presenti). cur denota il dominio corrente, next denota il dominio
 * che conterra' il nuovo valore delle energie. Questa funzione
 * realizza il passo 2 descritto nella specifica del progetto.
 * -------------------------------------------------------------
 * Modifica apportata:
 * per calcolare quali celle avessero il valore superiore a EMAX ho scelto di utilizzare la riduzione di OpenMP,
 * per gestire la propagazione e ogni cella conrolla quelle vicino a se per fare questo ho usato il pattern Stencil,
 * che ha comportato la creazione di una ghost area attorno alla matrice iniziale, questo ha evitato di fare ulteriori
 * controlli.
 */
void propagate_energy(float *cur, float *next, int n) {
  const float FDELTA = EMAX / 4;
  #pragma omp parallel for collapse(2) default (none) shared(cur, next, n)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      float F = *IDX(cur, i, j, n);
      float *out = IDX(next, i, j, n);

      /* Se l'energia del vicino di sinistra (se esiste) e'
         maggiore di EMAX, allora la cella (i,j) ricevera'
         energia addizionale FDELTA = EMAX/4 */
      if ( *IDX(cur, i, j - 1, n) > EMAX) {
        F += FDELTA;
      }
      /* Idem per il vicino di destra */
      if ( *IDX(cur, i, j + 1, n) > EMAX) {
        F += FDELTA;
      }
      /* Idem per il vicino in alto */
      if ( *IDX(cur, i - 1, j, n) > EMAX) {
        F += FDELTA;
      }
      /* Idem per il vicino in basso */
      if ( *IDX(cur, i + 1, j, n) > EMAX) {
        F += FDELTA;
      }

      if (F > EMAX) {
        F -= EMAX;
      }

      /* Si noti che il valore di F potrebbe essere ancora
         maggiore di EMAX; questo non e' un problema:
         l'eventuale eccesso verra' rilasciato al termine delle
         successive iterazioni vino a riportare il valore
         dell'energia sotto la foglia EMAX. */
      *out = F;
    }
  }
}

/**
 * Restituisce l'energia media delle celle del dominio grid di
 * dimensioni n*n. Il dominio non viene modificato.
 * ------------------------------------------------------------
 * Modificha apportata:
 * parallelizzazione con collasso (2) e aggiunta della clausola reduction 
 */
float average_energy(float *grid, int n) {
  float sum = 0.0f;
  #pragma omp parallel for collapse(2) default (none) shared(grid, n) reduction(+: sum)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      sum += *IDX(grid, i, j, n);
    }
  }
  return (sum / (n * n));
}

int main(int argc, char *argv[]) {
  float *cur, *next;
  int s, n = 256, nsteps = 2048;
  int m = n + 2 * HALO;
  //int num_threads = omp_get_num_threads();
  //FILE *omp_res = fopen("omp_out_"+argv[3], "a");

  srand(19); /* Inizializzazione del generatore pseudocasuale */

  if (argc > 4) {
    fprintf(stderr, "Usage: %s [nsteps] [n] [num_threads_used]\n", argv[0]);
    return EXIT_FAILURE;
  }

  if (argc > 1) {
    nsteps = atoi(argv[1]);
  }

  if (argc > 2) {
    n = atoi(argv[2]);
  }

  const size_t size = m * m * sizeof(float);

  /* Allochiamo i domini */
  cur = (float*) malloc(size);
  assert(cur);
  next = (float*) malloc(size);
  assert(next);

  /* L'energia iniziale di ciascuna cella e' scelta 
     con probabilita' uniforme nell'intervallo [0, EMAX*0.1] */
  setup(cur, n, 0, EMAX * 0.1);

  const double tstart = hpc_gettime();
  for (s = 0; s < nsteps; s++) {
    /* L'ordine delle istruzioni che seguono e' importante */
    increment_energy(cur, n, EDELTA);
    count_cells(cur, n);
    propagate_energy(cur, next, n);
    average_energy(next, n);
    float * tmp = cur;
    cur = next;
    next = tmp;
  }
  const double elapsed = hpc_gettime() - tstart;

  double Mupdates = (((double) n) * n / 1.0e6) * nsteps; /* milioni di celle aggiornate per ogni secondo di wall clock time */
  fprintf(stderr, "%s : %.4f Mupdates in %.4f seconds (%f Mupd/sec)\n", argv[0], Mupdates, elapsed, Mupdates / elapsed);
  //fprintf(omp_result,"%.4f\n", elapsed);
  //fclose(omp_result);
  //fprintf(stderr, "%.4f\n",elapsed);
  if (argc > 3) {
    FILE * omp_result = fopen(argv[3], "a");
    fprintf(omp_result, "%.4f\n", elapsed);
    fclose(omp_result);
  }

  /* Libera la memoria */
  free(cur);
  free(next);

  return EXIT_SUCCESS;
}
