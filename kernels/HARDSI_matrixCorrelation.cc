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
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the matrix correlation.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.
 */


// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int N;
MATRIX_DEFINE(DATA_TYPE, data);
MATRIX_DEFINE(DATA_TYPE, corr);

DATA_TYPE eps = 0.1;


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

	MATRIX_ALLOCATE(DATA_TYPE, N, N, data,		initialValue);
	MATRIX_ALLOCATE(DATA_TYPE, N, N, corr,		initialValue);
}


/**
 * \brief Matrix correlation
 */
void Worker::compute ()
{
	unsigned int i, j, k;
	float float_n = (float)N;
	DATA_TYPE mean[N], stddev[N];

	for (j = 0; j < N; j++)
	{
		mean[j] = 0.0;
		for (i = 0; i < N; i++)
		{
			mean[j] += MATRIX_GET(data, i, j);
		}
		mean[j] /= (float)N;
	}

	for (j = 0; j < N; j++)
	{
		stddev[j] = 0.0;
		for (i = 0; i < N; i++)
		{
			DATA_TYPE tmp = MATRIX_GET(data, i, j);
			stddev[j] += (tmp - mean[j]) * (tmp - mean[j]);
		}

		stddev[j] /= float_n;
		stddev[j] = sqrt(stddev[j]);

		// The following in an inelegant but usual way to handle near-zero std. dev. values, which below would cause a zero-divide.
		stddev[j] = stddev[j] <= eps ? 1.0 : stddev[j];
	}

	// Center and reduce the column vectors.
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			DATA_TYPE tmp = MATRIX_GET(data, i, j);
			tmp -= mean[j];
			tmp /= sqrt(float_n) * stddev[j];
			MATRIX_SET(data, i, j, tmp);
		}
	}

	// Calculate the m * m correlation matrix.
	for (i = 0; i < N; i++)
	{
		MATRIX_SET(corr, i, i, 1.0);
		for (j = i+1; j < N; j++)
		{
			MATRIX_SET(corr, i, j, 0.0);
			for (k = 0; k < N; k++)
			{
				DATA_TYPE tmp_i = MATRIX_GET(data, k, i);
				DATA_TYPE tmp_j = MATRIX_GET(data, k, j);
				MATRIX_ADD(corr, i, j,  tmp_i * tmp_j);
			}
//			corr[j][i] = corr[i][j];
		}
	}
	MATRIX_SET(corr, N-1, N-1, 1.0);
}


void Worker::postProcess()
{
	MATRIX_FREE(data,		N, N, DATA_TYPE);
	MATRIX_FREE(corr,		N, N, DATA_TYPE);
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
