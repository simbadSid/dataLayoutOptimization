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
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the matrix multiplication.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.
 */


// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int srcMatrix_X0, srcMatrix_Y0;
MATRIX_DEFINE(DATA_TYPE, srcMatrix0);

unsigned int srcMatrix_X1, srcMatrix_Y1;
MATRIX_DEFINE(DATA_TYPE, srcMatrix1);

MATRIX_DEFINE(DATA_TYPE, dstMatrix);


// --------------------------------------
// Local functions
// --------------------------------------
void Worker::preProcess(char *fileInput_name, char *fileOutput_name)
{
	InputData_matrix			*inputData	= new InputData_matrix(fileInput_name);
	InputData_matrix_Attribute	*attribute	= NULL;

	inputData->parseInputDataFile_matrixHeader(&srcMatrix_X0, &srcMatrix_Y0);
	parseInputDataFile_matrixBody(srcMatrix0, inputData, srcMatrix_X0, srcMatrix_Y0, attribute);	// Allocates the matrix 0 and fill it with data

	inputData->parseInputDataFile_matrixHeader(&srcMatrix_X1, &srcMatrix_Y1);
	parseInputDataFile_matrixBody(srcMatrix1, inputData, srcMatrix_X1, srcMatrix_Y1, attribute);	// Allocates the matrix 1 and fill it with data

	if (srcMatrix_X0 != srcMatrix_Y1)
	{
		LOGGER_error(1, "The parsed matrixes can't be multiplied: Matrix(%u, %u) * Matrix(%u, %u)", srcMatrix_X0, srcMatrix_Y0, srcMatrix_X1, srcMatrix_Y1);
	}

	MATRIX_ALLOCATE(DATA_TYPE, srcMatrix_X1, srcMatrix_Y0, dstMatrix, 0);							// Allocates dst matrix and fill it with 0

	delete (inputData);
}


/**
 * \brief Matrix multiplication based on the naive algorithm
 */
void Worker::compute ()
{
	unsigned int x, y, k;
	DATA_TYPE a0, a1, v;

	for (y=0; y<srcMatrix_Y0; ++y)
	{
		for (x=0; x<srcMatrix_X1; ++x)
		{
			v = 0;
			for (k=0; k<srcMatrix_X0; ++k)
			{
				a0 = MATRIX_GET(srcMatrix0, k, y);
				a1 = MATRIX_GET(srcMatrix1, x, k);
				v	+= a0*a1;
			}
			MATRIX_SET(dstMatrix, x, y, v);
		}
	}
}


void Worker::postProcess()
{
	unsigned int x, y;
	DATA_TYPE v;

	LOGGER(LOGGER_FLAG_INPUT_DATA, "Result of the multiplication\n");

	for (y=0; y<srcMatrix_Y0; ++y)
	{
		for (x=0; x<srcMatrix_X1; ++x)
		{
			v = MATRIX_GET(dstMatrix, x, y);

			// Done to avoid LOGGER -Werror
			v++;
			v--;
			LOGGER(LOGGER_FLAG_INPUT_DATA, " %f", (float)v);
		}
		LOGGER(LOGGER_FLAG_INPUT_DATA, "\n");
	}
	MATRIX_FREE(srcMatrix0,	srcMatrix_X0,	srcMatrix_Y0, DATA_TYPE);
	MATRIX_FREE(srcMatrix1,	srcMatrix_X1,	srcMatrix_Y1, DATA_TYPE);
	MATRIX_FREE(dstMatrix,	srcMatrix_X1,	srcMatrix_Y0, DATA_TYPE);
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
