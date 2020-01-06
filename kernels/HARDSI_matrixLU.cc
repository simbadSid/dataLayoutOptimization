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
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the matrix LU algorithm.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.
 */


// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int N;
MATRIX_DEFINE(DATA_TYPE, A);
MATRIX_DEFINE(DATA_TYPE, A_l);
MATRIX_DEFINE(DATA_TYPE, A_c);



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

	MATRIX_ALLOCATE(DATA_TYPE, N, N, A,		initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, A_l,	initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, A_c,	initialValue);
}


/**
 * \brief Matrix correlation
 */
void Worker::compute ()
{
	unsigned int i, j;
	unsigned int k;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j <i; j++)
		{
			DATA_TYPE a_ij = MATRIX_GET(A, i, j);
			for (k = 0; k < j; k++)
			{
				DATA_TYPE tmp_ik = MATRIX_GET(A_c,	i, k);
				DATA_TYPE tmp_kj = MATRIX_GET(A_l,	k, j);

				a_ij -= tmp_ik * tmp_kj;
			}
			a_ij /= MATRIX_GET(A_l,	j, j);
			MATRIX_SET(A, i, j, a_ij);
		}

		for (j = i; j < N; j++)
		{
			DATA_TYPE a_ij = MATRIX_GET(A, i, j);
			for (k = 0; k < i; k++)
			{
				DATA_TYPE tmp_ik = MATRIX_GET(A_c,	i, k);
				DATA_TYPE tmp_kj = MATRIX_GET(A_l,	k, j);
				a_ij -= tmp_ik * tmp_kj;
			}
			MATRIX_SET(A, i, j, a_ij);
		}
	}
}


void Worker::postProcess()
{
	MATRIX_FREE(A,		N, N, DATA_TYPE);
	MATRIX_FREE(A_l,	N, N, DATA_TYPE);
	MATRIX_FREE(A_c,	N, N, DATA_TYPE);
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

