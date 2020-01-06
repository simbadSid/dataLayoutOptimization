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
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the floyd warshall algorithm.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.
 */


// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int N;
MATRIX_DEFINE(DATA_TYPE, path);
MATRIX_DEFINE(DATA_TYPE, path_l);
MATRIX_DEFINE(DATA_TYPE, path_c);



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

	MATRIX_ALLOCATE(DATA_TYPE, N, N, path,		initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, path_l,	initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, path_c,	initialValue);
}


/**
 * \brief Matrix correlation
 */
void Worker::compute ()
{
	unsigned int i, j;
	unsigned int k;

	for (k = 0; k < N; k++)
//	for (k = 0; k < 1; k++)
	{
		for (j = 0; j < N; j++)
		{
			DATA_TYPE path_kj = MATRIX_GET(path_c,	k, j);
			for (i = 0; i < N; i++)
			{
				DATA_TYPE path_ik = MATRIX_GET(path_l, i, k);
				DATA_TYPE path_ij = MATRIX_GET(path,	i, j);

				MATRIX_SET(path, i, j, (path_ij < path_ik + path_kj)? path_ij : path_ik + path_kj);
			}
		}
	}
}


void Worker::postProcess()
{
	MATRIX_FREE(path,	N, N, DATA_TYPE);
	MATRIX_FREE(path_l,	N, N, DATA_TYPE);
	MATRIX_FREE(path_c,	N, N, DATA_TYPE);
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

