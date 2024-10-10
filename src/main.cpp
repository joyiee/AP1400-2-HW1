
#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"

int main(int argc, char **argv)
{
    if (false) // make false to run unit-tests
    {
        // debug section
        Matrix matrix = {{-1, 1.5, -1.75, -2}, {-2, 2.5, -2.75, -3}, {3, 3.5, -3.75, -4}, {4, 4.5, 4.75, -5}};
        //Matrix matrix = {{-12, 1.5, -1.75, -2}, {-2, 2.5, -2.75, -3}, {25, 3.5, -3.75, -4}, {4, 4.5, 4.75, -51}};
        Matrix inverse = {algebra::inverse(matrix)};
        algebra::show(inverse);
    }
    else
    {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}