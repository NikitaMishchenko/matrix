#include <iostream>

#include "matrix.h"

#include <string>


int main()
{
    matrix m;
    m.Load_matrix("A.txt");
    std::cout << m << std::endl;
    matrix m2 = m;
    m2.erase_column(0);
    m2.erase_row(1);
    m2.info();
    std::cout << m2 << std::endl;

    m.replace_data_part(m2,1,2,1,2);
    std::cout << m << std::endl;



}
