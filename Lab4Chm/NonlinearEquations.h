#pragma once
#include <fstream>
#include <string>
#include "../../../2/Lab2Chm/Lab2Chm/Matrix.h"
#include "../../../2/Lab2Chm/Lab2Chm/Vector.h"

#include "PolStr.h"

namespace luMath
{
	template<class T>
	class NonlinearEquations
	{
    public:
        enum class METHOD
        {
            NEWTON=1,
            ITERATIONS,
            STEEPESTDESCENT,
        };
    private:
	    METHOD              _method; // Требуемый метод решения СНУ
        int                 n;       // Размерность СНУ
        Vector<T>           x0;      // Начальное приближение
        T                   eps;     // Требуемая погрешность решения
        Vector<std::string> FunSys;  // Система функций
        std::ofstream       _fout;   // Выходной файл
    public:

        NonlinearEquations()
        {
            std::ifstream _fin("input.txt");
            _fout = std::ofstream("output.txt");
            int c;
            _fin >> c >> n;
            _method = static_cast<METHOD>(c);
            x0 = Vector<T>(n);
            x0.transposition();
            for (int i = 0; i < n; i++)
            {
                _fin >> c;
                x0[i] = c;
            }
            _fin >> eps;
            FunSys = Vector<std::string>(n);
            
            _fin.seekg(2, std::ios_base::cur);
            for (int i = 0; i < n; i++)
                getline(_fin, FunSys[i]);
                
            FunSys.transposition();
            std::cout << FunSys << '\n';
            _fin.close();
        }
        ~NonlinearEquations()
        {
            _fout.close();
        }

        METHOD getMethod() { return _method; }

        static Vector<T> GaussMethod(const Matrix<T>& A, const Vector<T>& b)
        {
            Matrix<T> expandedMatrix(A.getRows(), A.getCols() + 1);
            for (int i = 0; i < expandedMatrix.getRows(); i++)
                for (int j = 0; j < expandedMatrix.getCols(); j++)
                    if (j == expandedMatrix.getCols() - 1)
                        expandedMatrix[i][j] = b[i];
                    else
                        expandedMatrix[i][j] = A[i][j];
            return GaussMethod(expandedMatrix);
        }

        static Vector<T> GaussMethod(const Matrix<T>& expandedMatrix)
        {
            Matrix<T> tempMatrix(expandedMatrix);
           
            for (int i = 0; i < tempMatrix.getRows(); i++) 
            {
                T coeff = tempMatrix[i][i]; 
                

                for (int j = i; j < tempMatrix.getRows() + 1; j++)
                    tempMatrix[i][j] /= coeff;
                
                for (int j = i + 1; j < tempMatrix.getRows(); j++)
                {
                    coeff = tempMatrix[j][i]; 
                    for (int k = i; k < tempMatrix.getCols(); k++) 
                        tempMatrix[j][k] -= coeff * tempMatrix[i][k]; 
                }
            }

            Vector<T> result(tempMatrix.getRows());
            result.transposition();
            
            for (int i = result.getLength() - 1; i >= 0; i--)
            {
                T sumCoeff = 0;
                for (int j = i + 1; j < result.getLength(); j++)
                    sumCoeff += tempMatrix[i][j] * result[j];

                result[i] = tempMatrix[i][result.getLength()] - sumCoeff;
            }
            return result;
        }

        static Matrix<T> getInverseMatrixByMethod(Vector<T>(*Method)(const Matrix<T>&, const Vector<T>&), Matrix<T> matrix)
        {
            int m = matrix.getCols();
            Matrix<T> inverseMatrix(m);
            Vector<Vector<T>> x_temp(m);
            Vector<Vector<T>> E(m);
            for (int i = 0; i < m; i++)
            {
                E[i] = Vector<T>(m);
                E[i][i] = 1;
                E[i].transposition();
            }

            for (int i = 0; i < m; i++)
                x_temp[i] = Method(matrix, E[i]);
                
            
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    inverseMatrix[i][j] = x_temp[j][i];
            std::cout << std::setw(10) << inverseMatrix;
            return inverseMatrix;
        }

        Matrix<T> getJacobi(Vector<T> x) 
        {
            Matrix<T> J(n);
            for (int i = 0; i < n; i++)
            {
                std::cout << FunSys[i].c_str() << "\n";
                const char* polStr = CreatePolStr(FunSys[i].c_str(), n);
                
                if (GetError() == ERR_OK)
                {
                    for (int j = 0; j < n; j++)
                    {
                        const double* x_p = x.getPointer();
                        J[i][j] = EvalPolStr(polStr, x_p, 1, j + 1);
                        std::cout << EvalPolStr(polStr, x_p, 1, j+1) << "\n";
                    }
                }
                else std::cerr << "Error: " << GetError();
            }

            std::cout << J;
            return J;
        }


        void Newton() 
        {
            Matrix<T> J = getJacobi(x0);
            Matrix<T> J_1 = getInverseMatrixByMethod(NonlinearEquations::GaussMethod, J);
        
        }
        void Iterations() 
        {
        
        
        }
        void SteepestDescent() 
        {
        
        
        
        }

	};
}

