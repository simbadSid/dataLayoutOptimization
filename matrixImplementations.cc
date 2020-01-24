

/**
 * \brief This module gathers a subset of the matrix implementations used by the HARDSI compiler.
 * \detail In this context, a data-structure implementation is a set of routines used to define, allocate access and free a data structure
 * \author Riyane Sid Lakhdar
 */





/*===========================================
 * \brief Matrix definition and accessers depending on the data layout
 ===========================================*/
	#define MATRIX_LINE_MAJ_DEFINE(type, name)						type **name
	#define MATRIX_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)		ALLOCATE_2D_ARRAY(type, Y, X, name, initVal);
	#define MATRIX_LINE_MAJ_GET(m, x, y)							m[y][x]
	#define MATRIX_LINE_MAJ_SET(m, x, y, val)						m[y][x]	= val
	#define MATRIX_LINE_MAJ_ADD(m, x, y, val)						m[y][x]	+= val
	#define MATRIX_LINE_MAJ_FREE(m, X, Y, type)						FREE_2D_ARRAY(m, Y, X, type)

	#define MATRIX_COLUMN_MAJ_DEFINE(type, name)					type **name
	#define MATRIX_COLUMN_MAJ_ALLOCATE(type, X, Y, name, initVal)	ALLOCATE_2D_ARRAY(type, X, Y, name, initVal);
	#define MATRIX_COLUMN_MAJ_GET(m, x, y)							m[x][y]
	#define MATRIX_COLUMN_MAJ_SET(m, x, y, val)						m[x][y]	= val
	#define MATRIX_COLUMN_MAJ_ADD(m, x, y, val)						m[x][y]	+= val
	#define MATRIX_COLUMN_MAJ_FREE(m, X, Y, type)					FREE_2D_ARRAY(m, X, Y, type)

	#define MATRIX_DIAG_MAJ_DEFINE(type, name)						type **name
	#define MATRIX_DIAG_MAJ_ALLOCATE(type, X, Y, name, initVal)		ALLOCATE_2D_ARRAY(type, 1+X+Y, 1+X+Y, name, initVal);	// The correct formula is sqrt(x*x + Y*Y).   But we use the triangular inequality to reduce the number of computations (and cause we have no sqrt func)
	#define DIAG_MAJ_X(x, y)										((x)==(y)) ? (x) : ((x) > (y)) ? (y)			: (x)
	#define DIAG_MAJ_Y(x, y)										((x)==(y)) ? (0) : ((x) > (y)) ? (1+x+x-y-y)	: (y+y-x-x)
	#define MATRIX_DIAG_MAJ_GET(m, x, y)							m[DIAG_MAJ_X(x, y)][DIAG_MAJ_Y(x, y)]
	#define MATRIX_DIAG_MAJ_SET(m, x, y, val)						m[DIAG_MAJ_X(x, y)][DIAG_MAJ_Y(x, y)] = val;
	#define MATRIX_DIAG_MAJ_ADD(m, x, y, val)						m[DIAG_MAJ_X(x, y)][DIAG_MAJ_Y(x, y)] += val
	#define MATRIX_DIAG_MAJ_FREE(m, X, Y, type)						FREE_2D_ARRAY(m, 1+X+Y, 1+X+Y, type)


/*===========================================
 * \brief Matrix (stencil) definition and accessors depending on the data layout
 ===========================================*/
	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)							MATRIX_LINE_MAJ_DEFINE(type, name)
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)			MATRIX_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y)								MATRIX_LINE_MAJ_GET(m, x, y)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, val)							MATRIX_LINE_MAJ_SET(m, x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, val)							MATRIX_LINE_MAJ_ADD(m, x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, type)							MATRIX_LINE_MAJ_FREE(m, X, Y, type)
	#define MATRIX_STENCIL_LINE_MAJ_MAP_REDUCE_SUM_3(m, x, y)					{													\
																					float tmp	= MATRIX_LINE_MAJ_GET(m, x-1,y)		\
																								+ MATRIX_LINE_MAJ_GET(m, x,	y)		\
																								+ MATRIX_LINE_MAJ_GET(m, x+1,y)		\
																								+ MATRIX_LINE_MAJ_GET(m, x,	y-1)	\
																								+ MATRIX_LINE_MAJ_GET(m, x,	y+1);	\
																					MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, tmp);		\
																				}


/*
#include <pmmintrin.h>
	template<class T>
	struct MATRIX_STENCIL
	{
		T		**matrixLine;
		T		**matrixColumn;
	};

	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)							struct MATRIX_STENCIL<type> *name
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)			name = (MATRIX_STENCIL<type>*)SAFE_MALLOC(sizeof(MATRIX_STENCIL<type>));	\
																				MATRIX_LINE_MAJ_ALLOCATE(type, X, Y, (name->matrixLine), initVal);			\
																				MATRIX_COLUMN_MAJ_ALLOCATE(type, X, Y, (name->matrixColumn), initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y)								MATRIX_LINE_MAJ_GET((m->matrixLine), x, y)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, val)							MATRIX_LINE_MAJ_SET((m->matrixLine), x, y, val);							\
																				MATRIX_COLUMN_MAJ_SET((m->matrixColumn), x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, val)							MATRIX_LINE_MAJ_ADD((m->matrixLine), x, y, val);							\
																				MATRIX_COLUMN_MAJ_ADD((m->matrixColumn), x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, type)							MATRIX_LINE_MAJ_FREE((m->matrixLine), X, Y, type);							\
																				MATRIX_COLUMN_MAJ_FREE((m->matrixColumn), X, Y, type);						\
																				free(m);
	#define MATRIX_STENCIL_LINE_MAJ_MAP_REDUCE_SUM_3(m, x, y)					{																			\
																					__m128	vectLine, vectColumn, vectTmp;									\
																					float *ptrMatrixLine 	= (float*)(m->matrixLine	[y] + x);			\
																					float *ptrMatrixColumn	= (float*)(m->matrixColumn	[x] + y);			\
																																							\
																					vectLine		= _mm_loadu_ps(ptrMatrixLine);							\
																					vectColumn		= _mm_loadu_ps(ptrMatrixColumn);						\
																					vectTmp			= _mm_add_ps(vectLine, vectColumn);						\
																																							\
																					vectTmp			= _mm_hadd_ps(vectTmp,	vectTmp);						\
																					vectTmp			= _mm_hadd_ps(vectTmp,	vectTmp);						\
																																							\
																					vectLine[1]		= vectTmp[3] - vectLine[1] - vectLine[3] - vectColumn[3];\
																					vectColumn[1]	= vectLine[1];											\
																																							\
																					_mm_storeu_ps(ptrMatrixLine,	vectLine);								\
																					_mm_storeu_ps(ptrMatrixColumn,	vectColumn);							\
																				}
*/
/*#
	template<class T>
	struct MATRIX_STENCIL
	{
		T **matrixLine;
		T **matrixColumn;
	};

	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)					struct MATRIX_STENCIL<type> *name
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)	name = (MATRIX_STENCIL<type>*)SAFE_MALLOC(sizeof(MATRIX_STENCIL<type>));	\
																		MATRIX_LINE_MAJ_ALLOCATE(type, X, Y, (name->matrixLine), initVal)			;\
																		MATRIX_COLUMN_MAJ_ALLOCATE(type, X, Y, (name->matrixColumn), initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y)						MATRIX_LINE_MAJ_GET((m->matrixLine), x, y)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, val)					MATRIX_LINE_MAJ_SET((m->matrixLine), x, y, val)								;\
																		MATRIX_COLUMN_MAJ_SET((m->matrixColumn), x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, val)					MATRIX_LINE_MAJ_ADD((m->matrixLine), x, y, val)								;\
																		MATRIX_COLUMN_MAJ_ADD((m->matrixColumn), x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, type)					MATRIX_LINE_MAJ_FREE((m->matrixLine), X, Y, type)							;\
																		MATRIX_COLUMN_MAJ_FREE((m->matrixColumn), X, Y, type)						;\
																		free(m);
	#define MATRIX_STENCIL_LINE_MAJ_MAP(m, x, y, xl, xc, xr, yt, yb)	xl = MATRIX_LINE_MAJ_GET(m->matrixLine,		x-1,	y)						;\
																		xc = MATRIX_LINE_MAJ_GET(m->matrixLine,		x,		y)						;\
																		xr = MATRIX_LINE_MAJ_GET(m->matrixLine,		x+1,	y)						;\
																		yt = MATRIX_COLUMN_MAJ_GET(m->matrixColumn,	x,		y-1)					;\
																		yb = MATRIX_COLUMN_MAJ_GET(m->matrixColumn,	x,		y+1)

*/

/*
	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)													MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_DEFINE(type, name)
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)				MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y, X_block, Y_block)										MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x, y, X_block, Y_block)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, X_block, Y_block, val)									MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_SET(m, x, y, X_block, Y_block, val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, X_block, Y_block, val)									MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_ADD(m, x, y, X_block, Y_block, val)
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, X_block, Y_block, type)								MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_FREE(m, X, Y, X_block, Y_block, type)
	#define MATRIX_STENCIL_LINE_MAJ_MAP(m, x, y,X_block, Y_block, xl, xc, xr, yt, yb)					xl = MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x-1,	y,	X_block, Y_block)	;\
																										xc = MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x,	y,	X_block, Y_block)	;\
																										xr = MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x+1,	y,	X_block, Y_block)	;\
																										yt = MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x,	y-1,X_block, Y_block)	;\
																										yb = MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x,	y+1,X_block, Y_block)
*/
/*
	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)													MATRIX_LINE_MAJ_DEFINE(type, name)
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)				MATRIX_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y, X_block, Y_block)										MATRIX_LINE_MAJ_GET(m, x, y)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, X_block, Y_block, val)									MATRIX_LINE_MAJ_SET(m, x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, X_block, Y_block, val)									MATRIX_LINE_MAJ_ADD(m, x, y, val)
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, X_block, Y_block, type)								MATRIX_LINE_MAJ_FREE(m, X, Y, type)
	#define MATRIX_STENCIL_LINE_MAJ_MAP(m, x, y,X_block, Y_block, xl, xc, xr, yt, yb)					MATRIX_LINE_MAJ_MAP(m, x, y,xl, xc, xr, yt, yb)
*/

/*
	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)							MATRIX_LINE_MAJ_DEFINE(type, name)
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)			MATRIX_LINE_MAJ_ALLOCATE(type, 3*X, Y, name, initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y)								MATRIX_LINE_MAJ_GET(m, 3*x, y)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, val)							MATRIX_LINE_MAJ_SET(m, 3*x,		y,		val);								\
																				MATRIX_LINE_MAJ_SET(m, 3*x+2,	y-1,	val);								\
																				MATRIX_LINE_MAJ_SET(m, 3*x+1,	y+1,	val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, val)							{																			\
																					float tmp = val + MATRIX_STENCIL_COLUMN_MAJ_GET(m, x, y);				\
																					MATRIX_LINE_MAJ_ADD(m, 3*x,		y,		tmp);							\
																					MATRIX_LINE_MAJ_ADD(m, 3*x+2,	y-1,	tmp);							\
																					MATRIX_LINE_MAJ_ADD(m, 3*x+1,	y+1,	tmp);							\
																				}
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, type)							MATRIX_COLUMN_MAJ_FREE(m,3*X, Y, type)
	#define MATRIX_STENCIL_LINE_MAJ_MAP_REDUCE_SUM_3(m, x, y)					{																			\
																					float tmp	= MATRIX_LINE_MAJ_GET(m, (3*(x-1)), y)						\
																								+ MATRIX_LINE_MAJ_GET(m, (3*x),		y)						\
																								+ MATRIX_LINE_MAJ_GET(m, (3*x + 1), y)						\
																								+ MATRIX_LINE_MAJ_GET(m, (3*x + 2), y)						\
																								+ MATRIX_LINE_MAJ_GET(m, (3*x + 3), y);						\
																					MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, tmp);								\
																				}
*/


/*
#define XXXL	6
#define XXXC	10


	#define MATRIX_STENCIL_LINE_MAJ_DEFINE(type, name)							MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_DEFINE(type, name)
	#define MATRIX_STENCIL_LINE_MAJ_ALLOCATE(type, X, Y, name, initVal)			MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_ALLOCATE(type, 3*X, Y, XXXL, XXXC, name, initVal)
	#define MATRIX_STENCIL_LINE_MAJ_GET(m, x, y)								MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, 3*x, y, XXXL, XXXC)
	#define MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, val)							MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_SET(m, 3*x, y, XXXL, XXXC,		val);			\
																				MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_SET(m, 3*x+2,	y-1, XXXL, XXXC,	val);			\
																				MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_SET(m, 3*x+1,	y+1, XXXL, XXXC,	val)
	#define MATRIX_STENCIL_LINE_MAJ_ADD(m, x, y, val)							{																			\
																					float tmp = val + MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, 3*x, y);			\
																					MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_ADD(m, 3*x,		y, XXXL, XXXC,		tmp);		\
																					MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_ADD(m, 3*x+2,	y-1, XXXL, XXXC,	tmp);			\
																					MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_ADD(m, 3*x+1,	y+1, XXXL, XXXC,	tmp);			\
																				}
	#define MATRIX_STENCIL_LINE_MAJ_FREE(m, X, Y, type)							MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_FREE(m,3*X, Y, XXXL, XXXC, type)
	#define MATRIX_STENCIL_LINE_MAJ_MAP_REDUCE_SUM_3(m, x, y)					{																			\
																					float tmp	= MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, (3*(x-1)), y, XXXL, XXXC)	\
																								+ MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, (3*x),		y, XXXL, XXXC)	\
																								+ MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, (3*x + 1), y, XXXL, XXXC)	\
																								+ MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, (3*x + 2), y, XXXL, XXXC)	\
																								+ MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, (3*x + 3), y, XXXL, XXXC);	\
																					MATRIX_STENCIL_LINE_MAJ_SET(m, x, y, tmp);								\
																				}
*/


/*===========================================
 * \brief Matrix (tilled) definition and accessors depending on the data layout.
 *===========================================*/
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_LINE_DEFINE(type, name)											type ***name
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_LINE_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)		ALLOCATE_2D_ARRAY_SMALL_BLOCK(type, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), (Y_block), (X_block), (Y)%(Y_block), (X)%(X_block) , name, (initVal))
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_LINE_GET(m, x, y, X_block, Y_block)								m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(y, x, Y_block, X_block)]
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_LINE_SET(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(y, x, Y_block, X_block)]	= val
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_LINE_ADD(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(y, x, Y_block, X_block)]	+= val
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_LINE_FREE(m, X, Y, X_block, Y_block, type)						FREE_2D_ARRAY_SMALL_BLOCK(m, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), type)

	#define MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_DEFINE(type, name)											type ****name
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)		ALLOCATE_2D_ARRAY_LARGE_BLOCK(type, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), (Y_block), (X_block), (Y)%(Y_block), (X)%(X_block) , name, (initVal))
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_GET(m, x, y, X_block, Y_block)								m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)]
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_SET(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)]	= val
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_ADD(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)]	+= val
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_LINE_FREE(m, X, Y, X_block, Y_block, type)						FREE_2D_ARRAY_LARGE_BLOCK(m, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), (Y_block), (Y)%(Y_block), type)

	#define MATRIX_LINE_MAJ_SMALL_BLOCK_COLUMN_DEFINE(type, name)										type ***name
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_COLUMN_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)	ALLOCATE_2D_ARRAY_SMALL_BLOCK(type, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), (X_block), (Y_block), (X)%(X_block), (Y)%(Y_block) , name, (initVal))
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_COLUMN_GET(m, x, y, X_block, Y_block)							m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)]
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_COLUMN_SET(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)]	= val
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_COLUMN_ADD(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)]	+= val
	#define MATRIX_LINE_MAJ_SMALL_BLOCK_COLUMN_FREE(m, X, Y, X_block, Y_block, type)					FREE_2D_ARRAY_SMALL_BLOCK(m, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), type)

	#define MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_DEFINE(type, name)										type ****name
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)	ALLOCATE_2D_ARRAY_LARGE_BLOCK(type, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), (X_block), (Y_block), (X)%(X_block), (Y)%(Y_block), name, (initVal))
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_GET(m, x, y, X_block, Y_block)							m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)]
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_SET(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)]	= val
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_ADD(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)]	+= val
	#define MATRIX_LINE_MAJ_LARGE_BLOCK_COLUMN_FREE(m, X, Y, X_block, Y_block, type)					FREE_2D_ARRAY_LARGE_BLOCK(m, (int)ceil((float)(Y)/(float)(Y_block)), (int)ceil((float)(X)/(float)(X_block)), (X)%(X_block), (X_block), type)

	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_COLUMN_DEFINE(type, name)										type ***name
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_COLUMN_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)	ALLOCATE_2D_ARRAY_SMALL_BLOCK(type, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), (X_block), (Y_block), (X)%(X_block), (Y)%(Y_block) , name, (initVal))
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_COLUMN_GET(m, x, y, X_block, Y_block)							m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)]
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_COLUMN_SET(m, x, y, X_block, Y_block, val)					m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)]	= val
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_COLUMN_ADD(m, x, y, X_block, Y_block, val)					m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)]	+= val
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_COLUMN_FREE(m, X, Y, X_block, Y_block, type)					FREE_2D_ARRAY_SMALL_BLOCK(m, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block))type)

	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_COLUMN_DEFINE(type, name)										type ****name
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_COLUMN_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)	ALLOCATE_2D_ARRAY_LARGE_BLOCK(type, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), (X_block), (Y_block), (X)%(X_block), (Y)%(Y_block) , name, (initVal))
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_COLUMN_GET(m, x, y, X_block, Y_block)							m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)]
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_COLUMN_SET(m, x, y, X_block, Y_block, val)					m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)]	= val
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_COLUMN_ADD(m, x, y, X_block, Y_block, val)					m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)]	+= val
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_COLUMN_FREE(m, X, Y, X_block, Y_block, type)					FREE_2D_ARRAY_LARGE_BLOCK(m, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), (X_block), (X)%(X_block), type)

	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_LINE_DEFINE(type, name)										type ***name
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_LINE_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)	ALLOCATE_2D_ARRAY_SMALL_BLOCK(type, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), (Y_block), (X_block), (Y)%(Y_block), (X)%(X_block) , name, (initVal))
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_LINE_GET(m, x, y, X_block, Y_block)							m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(y, x, Y_block, X_block)]
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_LINE_SET(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(y, x, Y_block, X_block)]	= val
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_LINE_ADD(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_SMALL_BLOCK_INDEX(y, x, Y_block, X_block)]	+= val
	#define MATRIX_COLUMN_MAJ_SMALL_BLOCK_LINE_FREE(m, X, Y, X_block, Y_block, type)					FREE_2D_ARRAY_SMALL_BLOCK(m, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), type)

	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_LINE_DEFINE(type, name)										type ****name
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_LINE_ALLOCATE(type, X, Y, X_block, Y_block, name, initVal)	ALLOCATE_2D_ARRAY_LARGE_BLOCK(type, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), (Y_block), (X_block), (Y)%(Y_block), (X)%(X_block) , name, (initVal))
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_LINE_GET(m, x, y, X_block, Y_block)							m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)]
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_LINE_SET(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)]	= val
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_LINE_ADD(m, x, y, X_block, Y_block, val)						m[BLOCK_LINE_INDEX(x, X_block)][BLOCK_LINE_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(y, Y_block)][BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)]	+= val
	#define MATRIX_COLUMN_MAJ_LARGE_BLOCK_LINE_FREE(m, X, Y, X_block, Y_block, type)					FREE_2D_ARRAY_LARGE_BLOCK(m, (int)ceil((float)(X)/(float)(X_block)), (int)ceil((float)(Y)/(float)(Y_block)), (Y_block), (Y)%(Y_block), type)


/*===========================================
 * Tools
 ===========================================*/
	#define ALLOCATE_2D_ARRAY(type, X, Y, name, initVal)	name = (type **)SAFE_MALLOC(X * sizeof(type*));								\
															for (int MY_COMPILER_V0_X=0; MY_COMPILER_V0_X<(int)X; MY_COMPILER_V0_X++)	\
															{																			\
																name[MY_COMPILER_V0_X] = (type *)SAFE_MALLOC(Y*sizeof(type));			\
																memset((void*) (name[MY_COMPILER_V0_X]), initVal, Y*sizeof(type));		\
															}


	#define FREE_2D_ARRAY(name, X, Y, type)					for (int MY_COMPILER_V0_X=0; MY_COMPILER_V0_X<(int)X; MY_COMPILER_V0_X++)	\
															{																			\
																free(name[MY_COMPILER_V0_X]);											\
															}																			\
															free(name);

/*===========================================
 * \brief Tools for matrix by block
 *===========================================*/
	#define BLOCK_LINE_INDEX(x, X_block)							((x) / (X_block))
	#define BLOCK_LINE_SMALL_BLOCK_INDEX(x, y, X_block, Y_block)	(((x)%(X_block)) * (Y_block)) + ((y)%(Y_block))
	#define BLOCK_LINE_LARGE_BLOCK_INDEX(x, X_block)				((x) % (X_block))

	#define BLOCK_SIZE_X(NB_X, X, X_BLOCK, X_BLOCK_LAST)			((X_BLOCK_LAST > 0) && (X == NB_X-1))	?	X_BLOCK_LAST	: X_BLOCK
	#define BLOCK_SIZE_Y(NB_Y, Y, Y_BLOCK, Y_BLOCK_LAST)			((Y_BLOCK_LAST > 0) && (Y == NB_Y-1))	?	Y_BLOCK_LAST	: Y_BLOCK

	/**
	 * @param NB_X, NB_Y: number of block (per line and per column)
	 * @param X_BLOCK, Y_BLOCK: number of element in each block (per line and per column)
	 * @param X_BLOCK_LAST, Y_BLOCK_LAST: number of element in a block of the last line/column of blocks
	 */
	#define ALLOCATE_2D_ARRAY_SMALL_BLOCK(type, NB_X, NB_Y, X_BLOCK, Y_BLOCK, X_BLOCK_LAST, Y_BLOCK_LAST, name, initVal)														\
			name = (type ***)SAFE_MALLOC(NB_X * sizeof(type**));																												\
			for (int MY_COMPILER_V0_X=0; MY_COMPILER_V0_X<NB_X; MY_COMPILER_V0_X++)																								\
			{																																									\
				name[MY_COMPILER_V0_X] = (type **)SAFE_MALLOC(NB_Y * sizeof(type*));																							\
				for (int MY_COMPILER_V0_Y=0; MY_COMPILER_V0_Y<NB_Y; MY_COMPILER_V0_Y++)																							\
				{																																								\
					int MY_COMPILER_V0_BLOCK_SIZE = (BLOCK_SIZE_X(NB_X, MY_COMPILER_V0_X, X_BLOCK, X_BLOCK_LAST)) * (BLOCK_SIZE_Y(NB_Y, MY_COMPILER_V0_Y, Y_BLOCK, Y_BLOCK_LAST));	\
					name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y] = (type *)SAFE_MALLOC(MY_COMPILER_V0_BLOCK_SIZE * sizeof(type));													\
					memset((void*) (name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y]), initVal, MY_COMPILER_V0_BLOCK_SIZE * sizeof(type));												\
				}																																								\
			}

	/**
	 * @param NB_X, NB_Y: number of block (per line and per column)
	 * @param X_BLOCK, Y_BLOCK: number of element in each block (per line and per column)
	 * @param X_BLOCK_LAST, Y_BLOCK_LAST: number of element in a block of the last line/column of blocks
	 */
	#define ALLOCATE_2D_ARRAY_LARGE_BLOCK(type, NB_X, NB_Y, X_BLOCK, Y_BLOCK, X_BLOCK_LAST, Y_BLOCK_LAST, name, initVal)									\
			name = (type ****)SAFE_MALLOC(NB_X * sizeof(type***));																							\
			for (int MY_COMPILER_V0_X=0; MY_COMPILER_V0_X<NB_X; MY_COMPILER_V0_X++)																			\
			{																																				\
				name[MY_COMPILER_V0_X] = (type ***)SAFE_MALLOC(NB_Y * sizeof(type**));																		\
				for (int MY_COMPILER_V0_Y=0; MY_COMPILER_V0_Y<NB_Y; MY_COMPILER_V0_Y++)																		\
				{																																			\
					int MY_COMPILER_V0_BLOCK_SIZE_X = BLOCK_SIZE_X(NB_X, MY_COMPILER_V0_X, X_BLOCK, X_BLOCK_LAST);											\
					name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y] = (type **)SAFE_MALLOC(sizeof(type*) * MY_COMPILER_V0_BLOCK_SIZE_X);							\
					for (int MY_COMPILER_V0_XB=0; MY_COMPILER_V0_XB < MY_COMPILER_V0_BLOCK_SIZE_X; MY_COMPILER_V0_XB++)										\
					{																																		\
						int MY_COMPILER_V0_BLOCK_SIZE_Y = BLOCK_SIZE_Y(NB_Y, MY_COMPILER_V0_Y, Y_BLOCK, Y_BLOCK_LAST);										\
						name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y][MY_COMPILER_V0_XB] = (type *)SAFE_MALLOC(sizeof(type) * MY_COMPILER_V0_BLOCK_SIZE_Y);		\
						memset((void*) (name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y][MY_COMPILER_V0_XB]), initVal, sizeof(type) * MY_COMPILER_V0_BLOCK_SIZE_Y);	\
					}																																		\
				}																																			\
			}

	/**
	 * @param NB_X, NB_Y: number of block (per line and per column)
	 */
	#define FREE_2D_ARRAY_SMALL_BLOCK(name, NB_X, NB_Y, type)																					\
			for (int MX_COMPILER_V0_X=0; MX_COMPILER_V0_X<NB_X; MX_COMPILER_V0_X++)																\
			{																																	\
				for (int MX_COMPILER_V0_Y=0; MX_COMPILER_V0_Y<NB_Y; MX_COMPILER_V0_Y++)															\
				{																																\
					free(name[MX_COMPILER_V0_X][MX_COMPILER_V0_Y]);																				\
				}																																\
				free(name[MX_COMPILER_V0_X]);																									\
			}																																	\
			free(name)

	/**
	 * @param NB_X, NB_Y: number of block (per line and per column)
	 * @param X_block, Y_block: size of each block
	 * @param X_BLOCK_LAST, Y_BLOCK_LAST: number of element in a block of the last line/column of blocks
	 */
	#define FREE_2D_ARRAY_LARGE_BLOCK(name, NB_X, NB_Y, X_BLOCK, X_BLOCK_LAST, type)															\
			for (int MY_COMPILER_V0_X=0; MY_COMPILER_V0_X<NB_X; MY_COMPILER_V0_X++)																\
			{																																	\
				for (int MY_COMPILER_V0_Y=0; MY_COMPILER_V0_Y<NB_Y; MY_COMPILER_V0_Y++)															\
				{																																\
					int MY_COMPILER_V0_BLOCK_SIZE_X = BLOCK_SIZE_X(NB_X, MY_COMPILER_V0_X, X_BLOCK, X_BLOCK_LAST);								\
					for (int MY_COMPILER_V0_XB=0; MY_COMPILER_V0_XB<MY_COMPILER_V0_BLOCK_SIZE_X; MY_COMPILER_V0_XB++)							\
					{																															\
						free(name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y][MY_COMPILER_V0_XB]);														\
					}																															\
					free(name[MY_COMPILER_V0_X][MY_COMPILER_V0_Y]);																				\
				}																																\
				free(name[MY_COMPILER_V0_X]);																									\
			}																																	\
			free(name)


	void *SAFE_MALLOC(size_t size)
	{
		void *res = malloc(size);
		if (res == NULL)
		{
			printf("safe_malloc: Error while allocating dynamic memory at line %d: malloc(%d)", (unsigned int)size, __LINE__);
			exit(0);
		}
		return res;
	}



