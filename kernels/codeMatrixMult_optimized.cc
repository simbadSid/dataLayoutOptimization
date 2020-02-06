void matrixMult()
{
  int **a;
  int **b;
  int **res;

  a = (int**)malloc(N1*sizeof(int*));
  for (int i=0; i<N1; i++)
    a[i]=(int*)malloc(N0*sizeof(int));

  b = (int**)malloc(N2*sizeof(int*));
  for (int i=0; i<N2; i++)
    b[i]=(int*)malloc(N0*sizeof(int));

  res = (int**)malloc(N1*sizeof(int*));
  for (int i=0; i<N1; i++)
    res[i]=(int*)malloc(N2*sizeof(int))

  for (int j=0; j<N1; j++)
  {
     for (int i=0; i<N2; i++)
     {
        for (int k=0; k<N0; k++)
        {
           int a = a[j][k];
           int b = b[i][k];
           res[j][i] += a*b;
        }
     }
  }

  for (i=0; i<N0; i++)  free(a[i]);
  free(a);

  for (i=0; i<N2; i++)  free(b[i]);
  free(b);

  for (i=0; i<N2; i++)  free(res[i]);
  free(res);
}
