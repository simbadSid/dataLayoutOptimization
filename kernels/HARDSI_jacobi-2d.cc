/*
 *
 *  Created on: 20 aout. 2019
 *      Author: Riyane Sid Lakhdar
 */


#include "userInputParser.h"
#include "worker.h"
#include "util.h"
#include <stdio.h>
#include <math.h>




/**
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the jacobi function.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.
 */
/*
 * language_v0_matrixStencil_simpleStride2_jacobi.cc: called jacobi-2d in the polybench benchmark
 *
 *  Created on: 20 aout. 2019
 *      Author: rs254702
 */


#include "userInputParser.h"
#include "worker.h"
#include "util.h"
#include <stdio.h>
#include <math.h>





/**
 * \brief This files gives an implementation of the unit test worker (C++ class Worker).
 * \detail Computes a matrix multiplication using my v0 language:
 *        - No data-layout implementation is specified for any of the matrices
 */


// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int N;
MATRIX_DEFINE(DATA_TYPE, A_l);
MATRIX_DEFINE(DATA_TYPE, A_c);
MATRIX_DEFINE(DATA_TYPE, B_l);
MATRIX_DEFINE(DATA_TYPE, B_c);



// --------------------------------------
// Local functions
// --------------------------------------
void Worker::preProcess(char *fileInput_name, char *fileOutput_name)
{
	int initialValue = 43;

	FILE *f = fopen(fileInput_name, "r");
	if (f == NULL)
	{
		LOGGER_error(1, "Can't open the input-data file \"%s\"", fileInput_name);
	}
	if (1 > fscanf(f, "%d", &N))
	{
		LOGGER_error(1, "Can't open the input-data file \"%s\"", fileInput_name);
	}

	MATRIX_ALLOCATE(DATA_TYPE, N, N, A_l,	initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, A_c,	initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, B_l,	initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, B_c,	initialValue);
}


/**
 * \brief Matrix correlation
 */
void Worker::compute ()
{
	unsigned int i, j;
	unsigned int t;

	for (t = 0; t < N; t++)
	{
		for (i = 1; i < N - 1; i++)
		{
			DATA_TYPE tmpij_1;
			DATA_TYPE tmpij		= MATRIX_GET(A_c, i,	0);
			DATA_TYPE tmpij1	= MATRIX_GET(A_c, i,	1);
			for (j = 1; j < N - 1; j++)
			{
				tmpij_1	= tmpij;
				tmpij	= tmpij1;
				tmpij1	= MATRIX_GET(A_c, i, j+1);

				DATA_TYPE tmpi_1j	= MATRIX_GET(A_l, i-1,	j);
				DATA_TYPE tmpi1j	= MATRIX_GET(A_l, i+1,	j);

				MATRIX_SET(B_l, i, j, 0.2 * (tmpij + tmpij_1 + tmpij1 + tmpi1j + tmpi_1j));
				MATRIX_SET(B_c, i, j, 0.2 * (tmpij + tmpij_1 + tmpij1 + tmpi1j + tmpi_1j));
			}
		}

		for (i = 1; i < N - 1; i++)
		{
			DATA_TYPE tmpij_1;
			DATA_TYPE tmpij		= MATRIX_GET(B_c, i,	0);
			DATA_TYPE tmpij1	= MATRIX_GET(B_c, i,	1);
			for (j = 1; j < N - 1; j++)
			{
				tmpij_1	= tmpij;
				tmpij	= tmpij1;
				tmpij1	= MATRIX_GET(B_c, i, j+1);

				DATA_TYPE tmpi_1j	= MATRIX_GET(B_l, i-1,	j);
				DATA_TYPE tmpi1j	= MATRIX_GET(B_l, i+1,	j);

				MATRIX_SET(A_l, i, j, 0.2 * (tmpij + tmpij_1 + tmpij1 + tmpi1j + tmpi_1j));
				MATRIX_SET(A_c, i, j, 0.2 * (tmpij + tmpij_1 + tmpij1 + tmpi1j + tmpi_1j));
			}
		}
	}
}


void Worker::postProcess()
{
	MATRIX_FREE(A_l, N, N, DATA_TYPE);
	MATRIX_FREE(A_c, N, N, DATA_TYPE);
	MATRIX_FREE(B_l, N, N, DATA_TYPE);
	MATRIX_FREE(B_c, N, N, DATA_TYPE);
}

#ifdef MATRIX_OPTIMIZATION_STEP
/**
 * \brief Main function used by the compiler to run test executions.
 */
int main(int argc, char **argv)
{
	char fileInput_name[128], fileOutput_name[128], loggerFlag[128], trackerName[128];
	int cpuToPinMainThread;
	extractParameter(argc, argv, fileInput_name, fileOutput_name, loggerFlag, trackerName, &cpuToPinMainThread);
	LOGGER_Init(loggerFlag);
	LOGGER_parameters(LOGGER_FLAG_MAIN, *argv, (char*)MACROS_VALUE_STRING(DATA_TYPE), fileInput_name, fileOutput_name, trackerName, cpuToPinMainThread);
	if (cpuToPinMainThread >= 0)
		pinCurrentThreadOnCPU((unsigned)cpuToPinMainThread);

	Worker *worker = new Worker();
	worker->preProcess(fileInput_name, fileOutput_name);
	LOGGER(LOGGER_FLAG_MAIN, ">>>>>>> Beginning\n");

	worker->compute();

	LOGGER(LOGGER_FLAG_MAIN, ">>>>>>> End\n");
	worker->postProcess();

	delete worker;
}
#endif

