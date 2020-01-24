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
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the adi.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.
 */



// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int N;
MATRIX_DEFINE(DATA_TYPE, matrix);




#define STENCIL_FUNCTION(a0, a1, a2, a3, a4)	((a0*a1*a2) + (a3*a4))


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

	MATRIX_ALLOCATE(DATA_TYPE, N, N, matrix,	initialValue);
}


/**
 * \brief Matrix adi
 */
void Worker::compute ()
{
	unsigned int i, j;

	for (j = 1; j < N-1; j++)
	{
		for (i = 1; i < N-1; i++)
		{
			MATRIX_MAP_REDUCE_SUM_3(matrix, i, j);
		}
	}
}


void Worker::postProcess()
{
	MATRIX_FREE(matrix, N, N, DATA_TYPE);
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

