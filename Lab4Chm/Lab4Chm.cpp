#include <iostream>
#include "NonlinearEquations.h"

using namespace luMath;

int main()
{
    setlocale(LC_ALL, "Rus");
    NonlinearEquations<double> data;
    switch (data.getMethod())
    {
    case NonlinearEquations<double>::METHOD::NEWTON:  // Метод Ньютона
        data.Newton();
        break;
    case NonlinearEquations<double>::METHOD::ITERATIONS: // Метод Итераций
        data.Iterations();
        break;
    case NonlinearEquations<double>::METHOD::STEEPESTDESCENT: // Метод наискорейшего спуска
        data.SteepestDescent();
        break;
    }
    return 0;
}