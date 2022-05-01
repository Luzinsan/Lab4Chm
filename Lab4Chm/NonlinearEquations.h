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
        
        }
        void Iterations() 
        {
        
        
        }
        void SteepestDescent() 
        {
        
        
        
        }

	};
}

