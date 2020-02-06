void matrixMult()
{
  MATRIX_DEFINE(int, a);
  MATRIX_DEFINE(int, b);
  MATRIX_DEFINE(int, res);

  MATRIX_ALLOCATE(int, N0, N1, a);



  MATRIX_ALLOCATE(int, N2, N0, b);



  MATRIX_ALLOCATE(int, N2, N1, res);



  for (int j=0; j<N1; j++)
  {
     for (int i=0; i<N2; i++)
     {
        for (int k=0; k<N0; k++)
        {
           int a = MATRIX_GET(a, k, j);
           int b = MATRIX_GET(b, i, k);
           MATRIX_ADD(res, i,j, a*b);
        }
     }
  }

  MATRIX_FREE(a,   N0, N1, int);


  MATRIX_FREE(b,   N2, N0, int);


  MATRIX_FREE(res, N2, N1, int);

}
