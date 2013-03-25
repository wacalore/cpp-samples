#include "Matrix.h"
#include <cmath>
#include <iostream>

Vector operator+(const Vector &v,const Vector &w)
{
   int d=v.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = v[i]+w[i];
   return res;
}

Vector operator-(const Vector &v,const Vector &w)
{
   int d=v.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = v[i]-w[i];
   return res;
}

Vector operator*(const Vector &v,const Vector &w)
{
   int d=v.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = v[i]*w[i];
   return res;
}

double operator^(const Vector &v,const Vector &w)
{
   int d=v.size();
   double res=0.0;
   for(int i=0;i<d;i++) res = res + v[i]*w[i];
   return res;
}

Vector operator/(const Vector &v, const Vector &w)
{
    int d=v.size();
    Vector res(d);
    for(int i=0;i<d;i++) res[i] = v[i]/w[i];
    return res;
}


Vector operator+(const double &a,const Vector &w)
{
   int d=w.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = a+w[i];
   return res;
}

Vector operator-(const double &a,const Vector &w)
{
   int d=w.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = a-w[i];
   return res;
}

Vector operator*(const double &a,const Vector &w)
{
   int d=w.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = a*w[i];
   return res;
}

Vector operator/(const Vector &v,const double &a)
{
    int d=v.size();
    Vector res(d);
    for(int i=0;i<d;i++) res[i] = v[i]/a;
    return res;
}

Vector operator*(const Matrix &C, const Vector &v)
{
   int d=v.size();
   Vector res(d);
   for(int i=0;i<d;i++)
   {
      res[i] = 0.0;
      for(int k=0;k<d;k++) res[i] = res[i] + C[i][k]*v[k];
   }
   return res;
}

Matrix operator^(const Matrix& M, const Matrix& A)
{
    int dim = A[0].size();
    Matrix res(dim);
    Matrix transpose(dim);
    Matrix final(A.size());
    for(int l=0; l<dim; l++) res[l].resize(M.size());
    for(int i=0; i<dim; i++) transpose[i].resize(A.size());
    for(int i=0; i<A.size(); i++) final[i].resize(dim);

    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<A.size(); j++)
        {
            transpose[i][j] = A[j][i];
        }
    }

    for(int i = 0; i < dim; i++)
    {
        res[i] = M*transpose[i];
    }

    for(int i=0; i<A.size(); i++)
    {
        for(int j=0; j<dim; j++)
        {
            final[i][j] = res[j][i];
        }
    }

    return final;
}

Vector exp(const Vector &v)
{
   int d=v.size();
   Vector res(d);
   for(int i=0;i<d;i++) res[i] = exp(v[i]);
   return res;
}

Vector log(const Vector &v)
{
    int d=v.size();
    Vector res(d);
    for(int i=0; i<d; i++) res[i] = log(v[i]);
    return res;
}



double Det(Matrix &M)
{
    int dim = M[0].size();

    if( dim == 1 ) return M[0][0];

    double det = 0.0;

    for(int i = 0; i < dim; i++ )
    {
        Matrix Cofactor(dim-1);
        for(int l=0;l<dim-1;l++) Cofactor[l].resize(dim-1);

        int col=0, row=0;

        for(int j = 0; j < dim; j++)
        {
            if(j != i)
            {
                row = 0;
                for(int k = 1; k < dim; k++)
                {
                    Cofactor[row][col] = M[k][j];
                    row++;
                }
            col++;
            }
        }

        det += (i%2==1?-1.0:1.0) * M[0][i] * Det(Cofactor);
    }

    return det;
}

Matrix Inverse(Matrix &M)
{
    int dim = M.size();
    Matrix CofactorMatrix(dim);
    Matrix InverseMatrix(dim);

    for(int l=0;l<dim;l++) CofactorMatrix[l].resize(dim);
    for(int l=0;l<dim;l++) InverseMatrix[l].resize(dim);

    for(int i = 0; i < dim; i++)
    {
        for(int m = 0; m < dim; m++)
        {

            Matrix Cofactor(dim-1);
            for(int l=0;l<dim-1;l++) Cofactor[l].resize(dim-1);

            int col=0, row=0;
            for(int j = 0; j < dim; j++)
            {
                if(j != i)
                {
                    row = 0;
                    for(int k = 0; k < dim; k++)
                    {

                        if (k != m)
                        {

                            Cofactor[row][col] = M[k][j];
                            row++;
                        }
                    }
                col++;
                }
            }

            CofactorMatrix[i][m] = ((i+m)%2==1?-1.0:1.0) * Det(Cofactor);
        }
    }


    for(int i=0; i< dim; i++)
    {
        InverseMatrix[i] = (1.0/Det(M))*CofactorMatrix[i];
    }

    return InverseMatrix;
}

void DisplayMatrix(const Matrix& M)
{
    int rows = M.size();
    int columns = M[0].size();
    for(int i = 0; i < rows; i++)
    {
        std::cout << std::endl;
        for(int j=0; j<columns; j++)
        {
            std::cout << "M[" << i << "][" << j << "] = " << M[i][j] << "   ";
        }
    }
}





