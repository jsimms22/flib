#include "build3.hpp"

int main()
{
    constexpr int row = 2;
    constexpr int column = 2;
    std::array<double,row*column> mat = array_builder<double,row*column>(1.0);

    array_printer(row, column, mat);
}