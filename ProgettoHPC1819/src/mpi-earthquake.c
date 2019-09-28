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
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>     /* rand() */
#include <assert.h>

/* energia massima */
#define EMAX 4.0f
/* energia da aggiungere ad ogni timestep */
#define EDELTA 1e-4
#define HALO 1

/**
 * Restituisce un puntatore all'elemento di coordinate (i,j) del
 * dominio grid con n colonne.
 */
static inline float *IDX(float *grid, int i, int j, int n)
{
    return (grid + (i + HALO) * n + j + HALO);
}

/**
 * Restituisce un numero reale pseudocasuale con probabilita' uniforme
 * nell'intervallo [a, b], con a < b.
 */
float randab( float a, float b )
{
    return a + (b-a)*(rand() / (float)RAND_MAX);
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
void setup( float* grid, int n, float fmin, float fmax )
{
    for ( int i=0; i<n; i++ ) {
        for ( int j=0; j<n; j++ ) {
            *IDX(grid, i, j, n) = randab(fmin, fmax);
        }
    }
}

/**
 * Somma delta a tutte le celle del dominio grid di dimensioni
 * n*m. Questa funzione realizza il passo 1 descritto nella specifica
 * del progetto.
 */
void increment_energy( float *grid, int n, int m, float delta )
{
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            *IDX(grid, i, j, n) += delta;
        }
    }
}

/**
 * Restituisce il numero di celle la cui energia e' strettamente
 * maggiore di EMAX.
 */
int count_cells( float *grid, int n, int m)
{
    int c = 0;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            if ( *IDX(grid, i, j, n) > EMAX ) { c++; }
        }
    }
    return c;
}

/** 
 * Distribuisce l'energia di ogni cella a quelle adiacenti (se
 * presenti). cur denota il dominio corrente, next denota il dominio
 * che conterra' il nuovo valore delle energie. Questa funzione
 * realizza il passo 2 descritto nella specifica del progetto.
 */
void propagate_energy( float *cur, float *next, int n, int m)
{
    const float FDELTA = EMAX/4;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            float F = *IDX(cur, i, j, n);
            float *out = IDX(next, i, j, n);

            /* Se l'energia del vicino di sinistra (se esiste) e'
               maggiore di EMAX, allora la cella (i,j) ricevera'
               energia addizionale FDELTA = EMAX/4 */
            if ((j>0) && (*IDX(cur, i, j-1, n) > EMAX)) { F += FDELTA; }
            /* Idem per il vicino di destra */
            if ((j<n-1) && (*IDX(cur, i, j+1, n) > EMAX)) { F += FDELTA; }
            /* Idem per il vicino in alto */
            if ((i>0) && (*IDX(cur, i-1, j, n) > EMAX)) { F += FDELTA; }
            /* Idem per il vicino in basso */
            if ((i<n-1) && (*IDX(cur, i+1, j, n) > EMAX)) { F += FDELTA; }

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
 * Restituisce la somma dell'energia delle celle del dominio grid di
 * dimensioni n*m. Il dominio non viene modificato.
 */
float sum_energy(float *grid, int n, int m)
{
    float sum = 0.0f;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            sum += *IDX(grid, i, j, n);
        }
    }
    return sum;
}

int main( int argc, char* argv[] )
{   	
    float *cur, *next,Emean;
    int s,counter, n = 256, nsteps = 2048;
    int rank, comm_size;
    int m = n + 2 * HALO;
	
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
    MPI_Status status;

    srand(19); /* Inizializzazione del generatore pseudocasuale */
    
    if ( argc > 4 ) {
        fprintf(stderr, "Usage: %s [nsteps] [n] [num_threads_used]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if ( argc > 1 ) {
        nsteps = atoi(argv[1]);
    }

    if ( argc > 2 ) {
        n = atoi(argv[2]);
    }

    
    const size_t size = m*m*sizeof(float);

    MPI_Datatype rowtype;
    MPI_Type_vector(1, n, m, MPI_FLOAT, &rowtype);
    MPI_Type_commit(&rowtype);

    int *sendcnts = (int *)malloc(comm_size * sizeof(int)); 
	assert(sendcnts);
	sendcnts[0] = n / comm_size;
	
	int *displ = (int *)malloc(comm_size * sizeof(int)); 
	assert(displ);    
    displ[0] = 0;
    

    for (int i = 1; i < comm_size; i++) {
        displ[i] = displ[i - 1] + sendcnts[i - 1];
        sendcnts[i] = (i+1) * n / comm_size - i * n / comm_size;
    }

    /* Allochiamo i domini */
    cur = (float*)malloc(size); 
	assert(cur);
    next = (float*)malloc(size); 
	assert(next);

    /* L'energia iniziale di ciascuna cella e' scelta
       con probabilita' uniforme nell'intervallo [0, EMAX*0.1] */
    if (rank == 0) {
        setup(cur, n, 0, EMAX*0.1);
    }

    const double tstart = hpc_gettime();
	
    MPI_Scatterv(
		IDX(cur, 0, 0, n), 
		sendcnts, 
		displ, 
		rowtype, 
		IDX(cur, 0, 0, n), 
		m, 
		rowtype, 
		0, 
		MPI_COMM_WORLD);
		
    for (s=0; s<nsteps; s++) {
        /* L'ordine delle istruzioni che seguono e' importante */
        increment_energy(cur, n, sendcnts[rank], EDELTA);
        counter = 0;
        int local_c = count_cells(cur, n, sendcnts[rank]);
        MPI_Reduce(
			&local_c, 
			&counter,
			1, 
			MPI_INT, 
			MPI_SUM, 
			0, 
			MPI_COMM_WORLD);

    int upDest = rank - 1;
	int downDest = rank + 1;
	
	if (upDest < 0) {
	    upDest = MPI_PROC_NULL;
	}
	if (downDest == comm_size) {
	    downDest = MPI_PROC_NULL;
	}
        int sender = rank - 1;
        MPI_Sendrecv(
            IDX(cur, sendcnts[sender] - HALO, 0, n),
            HALO,
            rowtype,
            downDest,
            0,
            &cur[HALO],
            HALO,
            rowtype,
            upDest,
            0,
            MPI_COMM_WORLD,
            &status
        );
        MPI_Sendrecv(
            IDX(cur, 0, 0, n),
            HALO,
            rowtype,
            upDest,
            0,
            IDX(cur, sendcnts[sender], 0, n),
            HALO,
            rowtype,
            downDest,
            0,
            MPI_COMM_WORLD,
            &status
        );
        propagate_energy(cur, next, n, sendcnts[rank]);
        float sum_e = sum_energy(next, n, sendcnts[rank]);
        MPI_Reduce(
			&sum_e, 
			&Emean,
			1, 
			MPI_FLOAT, 
			MPI_SUM, 
			0, 
			MPI_COMM_WORLD);

        if (rank == 0) {
            Emean = (Emean / (n*n));
        }
        float *tmp = cur;
        cur = next;
        next = tmp;
    }
    const double elapsed = hpc_gettime() - tstart;
   
    if(argc > 3){
	if(rank==0){
		double Mupdates = (((double)n)*n/1.0e6)*nsteps; /* milioni di celle aggiornate per ogni secondo di wall clock time */
    		fprintf(stderr, "%s : %.4f Mupdates in %.4f seconds (%f Mupd/sec)\n", argv[0], Mupdates, elapsed, Mupdates/elapsed);
		FILE *mpi_out = fopen(argv[3],"a");
		fprintf(mpi_out,"%.4f\n",elapsed);
		fclose(mpi_out);
	}
    }

    /* Libera la memoria */
    free(cur);
    free(next);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
