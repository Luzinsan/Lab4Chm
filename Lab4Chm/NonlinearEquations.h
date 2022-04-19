#pragma once
#include <fstream>
#include <string>
#include "../../../2/Lab2Chm/Lab2Chm/Vector.h"

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
        METHOD _method;
        int                 n;      // ����������� ���
        T                   x0;     // ��������� �����������
        T                   eps;    // ��������� ����������� �������
        Vector<char*>*       FunSys; // ������� �������
        std::ofstream      _fout;   // �������� ����
    public:

        NonlinearEquations() 
        {
            std::ifstream _fin("input.txt");
            _fout = std::ofstream("output.txt");
            int method;
            _fin >> method >> n >> x0 >> eps;
            _method = static_cast<METHOD>(method);
           
            const char* str = "111";
            FunSys = new Vector<char*>(n, str, true);
           
            
            //for (int i = 0; i < n; i++)
            //   FunSys[i] = " ";
                
        
            std::cout << FunSys << '\n';
            _fin.close();
        }
        ~NonlinearEquations()
        {
            delete FunSys;
            _fout.close();
        }

        METHOD getMethod() { return _method; }

        void Newton() 
        {
        
        
        }
        void Iterations() 
        {
        
        
        }
        void SteepestDescent() 
        {
        
        
        
        }

	};
}

