/*
 *
 *  Created on: 25 févr. 2019
 *      Author: Riyane Sid Lakhdar
 */

#include <cmath>

#include "worker.h"
#include "inputData_jpeg.h"
#include "userInputParser.h"




/**
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the Jpeg Compression.
 * \detail This code is based on an original implementation of the JPEG compression: https://github.com/Ipsedo/CompressionJPEG.
 */


// --------------------------------------
// Local parameters (static)
// --------------------------------------
unsigned int srcImage_X, srcImage_Y;
unsigned int nbrBlock;


// --------------------------------------
// Pointers for each type of access on memory
// --------------------------------------
int 	**	srcImage;
double	***	srcImageBlock;
MATRIX_DEFINE(double,	srcBlock);
MATRIX_DEFINE(rle,		srcImageBlockPerLine);

//int		**	dstImage;
std::string dstImage;


// --------------------------------------
// Auxiliary functions (header)
// --------------------------------------
void		splitImg		();
double**	makeBlock		(double **srcBlock, unsigned int x, unsigned int y);
void		dctTransform	(double **srcBlock);
double		dct				(double **srcBlock, unsigned int x, unsigned int y);
void		quantify		(double **srcBlock);
void		zigZagEncoding	(double **srcBlock, unsigned int bl);
void 		rle_encode		(unsigned int bl);
void		dc_ac_encode	(unsigned int bl);
std::string	encode_huffman	(unsigned int bl);


#define GET_NBR_BLOCK(X, Y)		DIVIDE_ROUND_SUP(X, JPEG_BLOCK_SIZE) * DIVIDE_ROUND_SUP(Y, JPEG_BLOCK_SIZE)
#define GET_BLOCK_INDEX(x, y)	(y/JPEG_BLOCK_SIZE)*(srcImage_X/JPEG_BLOCK_SIZE) + (x/JPEG_BLOCK_SIZE)


// --------------------------------------
// Local functions
// --------------------------------------
void Worker::preProcess(char *fileInput_name, char *fileOutput_name)
{
	unsigned int i;
	InputData_jpeg				*inputData	= new InputData_jpeg(fileInput_name, fileOutput_name);
	InputData_jpeg_Attribute	*attribute	= NULL;

	inputData->parseInputDataFile_jpegHeader(&srcImage_X, &srcImage_Y);
	parseInputDataFile_jpegBody(srcImage, inputData, srcImage_X, srcImage_Y, attribute);			// Allocates the source matrix bitmap and fill it with data
	delete (inputData);

	nbrBlock = GET_NBR_BLOCK(srcImage_X, srcImage_Y);

//	MATRIX_ALLOCATE(int, srcImage_X, srcImage_Y, dstImage, 0);										// Allocates the destination matrix bitmap
	srcImageBlock = (double***)safe_malloc(sizeof(double**) * nbrBlock);							// Allocates a tilled version of the source image

	for (i=0; i<nbrBlock; i++)
	{
			MATRIX_ALLOCATE(double, JPEG_BLOCK_SIZE, JPEG_BLOCK_SIZE, srcBlock, 0);					// Allocates a matrix for each till of the input image
			srcImageBlock[i] = srcBlock;
	}
	MATRIX_ALLOCATE(struct RLE, nbrBlock, JPEG_BLOCK_NBR_ELEM, srcImageBlockPerLine, 0);					// Allocates the matrix containing one block (as a 1D array) per line
}


/**
 * \brief JPEG-compression of bitmap file
 */
void Worker::compute ()
{
	// Split the picture into blocks (and convert it to double)
	splitImg();

	unsigned int bl;

	// Apply the DCT transform to each block
	for (bl = 0; bl<nbrBlock; bl++)
	{
		dctTransform(srcImageBlock[bl]);
	}

	// Apply the quantification
	for (bl = 0; bl<nbrBlock; bl++)
	{
		quantify(srcImageBlock[bl]);
	}

	// Apply the zig-zag encoding and fill the dstMatrx (not tilled anymore)
	for (bl = 0; bl<nbrBlock; bl++)
	{
		zigZagEncoding(srcImageBlock[bl], bl);
	}

	// Apply the RLE compression
	for (bl=0; bl<nbrBlock; bl++)
	{
		rle_encode(bl);
	}

	// Code the coefficients
	for (bl = 0; bl < nbrBlock; bl++)
	{
		dc_ac_encode(bl);
	}

	// Codage de Huffman
	dstImage = "";
	for (bl = 0; bl < nbrBlock; bl++)
	{
		dstImage += encode_huffman(bl); //, tools.DC_code, tools.AC_code);
	}

}


void Worker::postProcess()
{
	unsigned int bl, i;
	rle dc;

	LOGGER(LOGGER_FLAG_INPUT_DATA, "Result of the multiplication\n");

	// TODO Write result image to file

	// Free the source image
	MATRIX_FREE(srcImage, srcImage_X, srcImage_Y, double);

	// Free the destinationimage
//	MATRIX_FREE(dstImage, srcImage_X, srcImage_Y, double);

	// Free each block of the tilled source image
	for (bl=0; bl<nbrBlock; bl++)
	{
		srcBlock = srcImageBlock[bl];
		MATRIX_FREE(srcBlock, JPEG_BLOCK_SIZE, JPEG_BLOCK_SIZE,	double);
	}
	// Free the whole block image
	free(srcImageBlock);

	// Free the matrix containing one block (as a 1D array) per line
	for (bl = 0; bl < nbrBlock; bl++)
	{
		for (i = 1; i < JPEG_BLOCK_NBR_ELEM; i++)
		{
			dc = MATRIX_GET(srcImageBlockPerLine, bl, i);
//TODO free(dc.binDcAc);
		}
	}
	MATRIX_FREE(srcImageBlockPerLine, nbrBlock, JPEG_BLOCK_NBR_ELEM, rle);
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


// --------------------------------------
// Auxiliary functions (header)
// --------------------------------------
void splitImg()
{
	unsigned int x, y;

	for (y=0; y<srcImage_Y; y+=JPEG_BLOCK_SIZE)
	{
		for (x=0; x<srcImage_X; x+=JPEG_BLOCK_SIZE)
		{
			makeBlock(srcImageBlock[GET_BLOCK_INDEX(x, y)], x, y);
		}
	}
}


double **makeBlock(double **srcBlock, unsigned int x, unsigned int y)
{
	int pixel;
	unsigned int xi, yi;
	unsigned int xPixel, yPixel;

	for (yi = 0; yi < JPEG_BLOCK_SIZE; yi++)
	{
		for (xi = 0; xi < JPEG_BLOCK_SIZE; xi++)
		{
			// Case within the picture
			if ((yi + y < srcImage_Y) && (xi + x < srcImage_X))
			{
				xPixel = x+xi;
				yPixel = y+yi;
			}
			// Case out of the picture: symmetrical padding
			else
			{
				xPixel = srcImage_X - xi - 1;
				yPixel = srcImage_Y - yi - 1;
			}
			pixel = MATRIX_GET(srcImage, xPixel, yPixel);

			MATRIX_SET(srcBlock, xi, yi, (double)pixel);
		}
	}

	return srcBlock;
}


void dctTransform(double **srcBlock)
{
	unsigned int xi, yi;
	double v;

	for (yi = 0; yi < JPEG_BLOCK_SIZE; yi++)
	{
		for (xi = 0; xi < JPEG_BLOCK_SIZE; xi++)
		{
			v = dct(srcBlock, xi, yi);
			MATRIX_SET(srcBlock, xi, yi, v);
		}
	}
}


double dct(double **srcBlock, unsigned int x, unsigned int y)
{
	double sum = 0., v;
	unsigned int xi, yi;

	for (xi = 0; xi < JPEG_BLOCK_SIZE; xi++)
	{
		for (yi = 0; yi < JPEG_BLOCK_SIZE; yi++)
		{
			v = MATRIX_GET(srcBlock, xi, yi);

			sum += v *	cos((2.0 * xi + 1.0) * x * M_PI / (2.0 * JPEG_BLOCK_SIZE)) *
						cos((2.0 * yi + 1.0) * y * M_PI / (2.0 * JPEG_BLOCK_SIZE));
		}
	}

	double c_i = (x == 0) ? 1.0 / sqrt(2.0) : 1.0;
	double c_j = (y == 0) ? 1.0 / sqrt(2.0) : 1.0;
	return 2 * c_i * c_j * sum / JPEG_BLOCK_SIZE;
}


void quantify(double **srcBlock)
{
	unsigned int xi, yi;
	double v;

	for (yi = 0; yi < JPEG_BLOCK_SIZE; yi++)
	{
		for (xi = 0; xi < JPEG_BLOCK_SIZE; xi++)
		{
			v = MATRIX_GET(srcBlock, xi, yi);
			v = round(v / QUANTIZATION_MATRIX[xi][yi]);
			MATRIX_SET(srcBlock, xi, yi, v);
		}
	}
}


void zigZagEncoding (double **srcBlock, unsigned int bl)
{
	unsigned xi, yi, nbSetElem, curr_diag;
	bool step;
	rle v;

	step		= true;											// True means up diagonal
	nbSetElem	= 0;
	xi			= 0;
	yi			= 0;
	curr_diag	= 0;

	while (nbSetElem < JPEG_BLOCK_NBR_ELEM)
	{
		for (;	xi >= 0 && xi < JPEG_BLOCK_SIZE &&
				yi >= 0 && yi < JPEG_BLOCK_SIZE;
				nbSetElem++)
		{
			v.nextNonNull		= 0;
			v.zeroAndMagnitude	= 0;
			v.pixel				= MATRIX_GET(srcBlock, xi, yi);
			MATRIX_SET(srcImageBlockPerLine, bl, nbSetElem, v);

			if (step)	{ xi++; yi--; }							// up right diagonal
			else		{ xi--; yi++; }							// down left diagonal
		}

		curr_diag++;

		if (curr_diag < JPEG_BLOCK_SIZE)						// Up triangle of the matrix
		{
			if (step)	yi++;
			else		xi++;
		}
		else													// Down triangle of the matrix
		{
			if (step)	{ xi--; yi += 2; }
			else		{ yi--; xi += 2; }
		}

		step = !step;
	}
}


void rle_encode(unsigned int bl)
{
	unsigned int i, j, nb_zero, zero_restant;
	bool all_zero;
	rle v, v0;

	// For each coeficient of the block
	for (i = 0; i < JPEG_BLOCK_NBR_ELEM; i++)
	{
		all_zero	= true;
		nb_zero		= 0;

		// Count the number of 0 (check if only zeros)
		for (j = i; j < JPEG_BLOCK_NBR_ELEM; j++)
		{
			v = MATRIX_GET(srcImageBlockPerLine, bl, j);
			if (abs(v.pixel) != 0)
			{
				all_zero = false;
				break;
			}
			else
			{
				nb_zero++;
			}
		}

		v = MATRIX_GET(srcImageBlockPerLine, bl, i);

		// Case only 0 in block: write EOB
		if (all_zero)
		{
			v.zeroAndMagnitude	= EOB;
			v.nextNonNull		= 0;
			MATRIX_SET(srcImageBlockPerLine, bl, i, v);
			break;
		}
		// Case not only 0 in block: add ZRL if more than 15 zeros
		else
		{
			zero_restant = nb_zero;
			if (zero_restant >= 16)
			{
				v.zeroAndMagnitude	= ZRL;
				v.nextNonNull		= 0;
				zero_restant		-= 16;
				i					+= 16;
			}

			v.zeroAndMagnitude = (unsigned char) zero_restant << 4;								// On indique le nombre de 0 sur les 4 premier bit de l'octet to_write
			i += zero_restant;																	// On incrémente i en fonction du nombre de zeros restant
			v0 = MATRIX_GET(srcImageBlockPerLine, bl, i);
			v.zeroAndMagnitude |= (unsigned char) ceil(log(abs(v0.pixel) + 1) / log(2.0));		// On calcule la magnitude du prochain coefficient et on l'ajoute au 4 derniers bits de l'octet to_write
			MATRIX_SET(srcImageBlockPerLine, bl, i, v);
		}
	}
}


void dc_ac_encode(unsigned int bl)
{
	unsigned int i;
	rle dc, ac;

	dc			= MATRIX_GET(srcImageBlockPerLine, bl, 0);
	dc.binDcAc	= (char*)malloc(sizeof(char) * dc.zeroAndMagnitude);
	write_bits<11>(dc.zeroAndMagnitude, dc.nextNonNull, dc.binDcAc);
	MATRIX_SET(srcImageBlockPerLine, bl, 0, dc);

	for (i = 1; i < JPEG_BLOCK_NBR_ELEM; i++)
	{
		ac = MATRIX_GET(srcImageBlockPerLine, bl,i);

		if (ac.zeroAndMagnitude == ZRL)
		{
			dc.pixel			= 0;
			dc.zeroAndMagnitude	= ZRL;
			dc.nextNonNull		= 0;
			dc.binDcAc			= (char*)malloc(sizeof(char) * 1);
			strcpy(dc.binDcAc, "\0");
			MATRIX_SET(srcImageBlockPerLine, bl, i, dc);
			continue;
		}

		if (ac.zeroAndMagnitude == EOB)
		{
			dc.pixel			= 0;
			dc.zeroAndMagnitude	= EOB;
			dc.nextNonNull		= 0;
			dc.binDcAc			= (char*)malloc(sizeof(char) * 1);
			strcpy(dc.binDcAc, "\0");
			MATRIX_SET(srcImageBlockPerLine, bl, i, dc);
			break;
		}

		dc.pixel			= 0;
		dc.zeroAndMagnitude	= ac.zeroAndMagnitude;
		dc.nextNonNull		= 0;
		dc.binDcAc			= (char*)malloc(sizeof(char) * (ac.zeroAndMagnitude & 0x0F));
		write_bits<10>(ac.zeroAndMagnitude & 0x0F, ac.nextNonNull, dc.binDcAc);
		MATRIX_SET(srcImageBlockPerLine, bl, i, dc);
	}
}


std::string encode_huffman(unsigned int bl)
{
// TODO free each dc.binDcAc
	rle dc, ac;
	unsigned int i, magn, zero_n_magn;
	std::string dc_bits_coeff_str, ac_bits_coeff_str;

	dc					= MATRIX_GET(srcImageBlockPerLine, bl, 0);
	magn				= dc.zeroAndMagnitude;
	dc_bits_coeff_str	= std::string(dc.binDcAc);

	if (magn == EOB)
	{
		return DC1_LENGTH[magn];
	}

	std::string res;
	res += DC1_LENGTH[magn] + dc_bits_coeff_str;

	for (i = 1; i < JPEG_BLOCK_NBR_ELEM; i++)
	{
		ac					= MATRIX_GET(srcImageBlockPerLine, bl, i);
		zero_n_magn			= ac.zeroAndMagnitude;
		ac_bits_coeff_str	= std::string(ac.binDcAc);

		if (zero_n_magn == EOB)
		{
			return res += AC_CODE[EOB];
		}

		if (zero_n_magn == ZRL)
		{
			res += AC_CODE[ZRL];
			continue;
		}

		res += AC_CODE[zero_n_magn] + ac_bits_coeff_str;
	}
	return res;
}

