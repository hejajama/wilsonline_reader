#include "unit_test_framework.hpp"

#include "ipglasma.hpp"
#include "wilsonline.hpp"

#include <iostream>
#include <vector>
#include <complex>
using namespace std;

TEST(MatrixProduct)
{
    std::vector<std::vector<std::complex<double>>> matrix1 =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };
    std::vector<std::vector<std::complex<double>>> matrix2 =  { {{2, 2}, {1, -1}, {0, 1}}, {{3, 0}, {2, 2}, {1, -1}}, {{1, 1}, {0, 0}, {2, -2}} };


    WilsonLine w1(matrix1);
    WilsonLine w2(matrix2);
    WilsonLine product = w1*w2;
    std::vector<std::vector<std::complex<double>>> correct_prodcut = {
        { {9,4}, {8,2}, {6,-8} },
        { {1,3}, {5,1}, {7,-2} },
        { {2,1}, {3,-3}, {3,0} }
    };

    // Test matrix product
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ASSERT_ALMOST_EQUAL(product.Element(i,j).real(), correct_prodcut[i][j].real(), 1e-10);
            ASSERT_ALMOST_EQUAL(product.Element(i,j).imag(), correct_prodcut[i][j].imag(), 1e-10);
        }
    }   
}

TEST(MultiplyByHermitianConjugate)
{
    std::vector<std::vector<std::complex<double>>> matrix1 =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };
    std::vector<std::vector<std::complex<double>>> matrix2 =  { {{2, 2}, {1, -1}, {0, 1}}, {{3, 0}, {2, 2}, {1, -1}}, {{1, 1}, {0, 0}, {2, -2}} };


    WilsonLine w1(matrix1);
    WilsonLine w2(matrix2);
    WilsonLine product = w1.MultiplyByHermitianConjugate(w2);
    std::vector<std::vector<std::complex<double>>> correct_prodcut = {
        { {7, -2}, {8, 0}, {8, 6} },
        { {6, 0}, {0, 3}, {1, 9} },
        { {4, -4}, {1, 0}, {1, 3} }
    };

    // Test matrix product
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ASSERT_ALMOST_EQUAL(product.Element(i,j).real(), correct_prodcut[i][j].real(), 1e-10);
            ASSERT_ALMOST_EQUAL(product.Element(i,j).imag(), correct_prodcut[i][j].imag(), 1e-10);
        }
    }   
}

TEST(WilsonLineAddition)
{
    std::vector<std::vector<std::complex<double>>> matrix1 =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };
    std::vector<std::vector<std::complex<double>>> matrix2 =  { {{2, 2}, {1, -1}, {0, 1}}, {{3, 0}, {2, 2}, {1, -1}}, {{1, 1}, {0, 0}, {2, -2}} };

    WilsonLine w1(matrix1);
    WilsonLine w2(matrix2);
    WilsonLine sum = w1 + w2;
    std::vector<std::vector<std::complex<double>>> correct_sum = {
        { {3, 3}, {3, -2}, {3, 1} },
        { {3, 1}, {3, 1}, {3, 1} },
        { {2, 1}, {0, -1}, {3, -1} }
    };

    // Test matrix addition
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ASSERT_ALMOST_EQUAL(sum.Element(i,j).real(), correct_sum[i][j].real(), 1e-10);
            ASSERT_ALMOST_EQUAL(sum.Element(i,j).imag(), correct_sum[i][j].imag(), 1e-10);
        }
    }   
}

TEST(WilsonLineTrace)
{
    std::vector<std::vector<std::complex<double>>> matrix =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };

    WilsonLine w(matrix);
    std::complex<double> trace = w.Trace();
    std::complex<double> correct_trace = {3, 1};

    // Test trace
    ASSERT_ALMOST_EQUAL(trace.real(), correct_trace.real(), 1e-10);
    ASSERT_ALMOST_EQUAL(trace.imag(), correct_trace.imag(), 1e-10);
}

TEST(WilsonLineComplexConjugate)
{
    std::vector<std::vector<std::complex<double>>> matrix =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };

    WilsonLine w(matrix);
    WilsonLine conjugate = w.ComplexConjugate();
    std::vector<std::vector<std::complex<double>>> correct_conjugate = {
        { {1, -1}, {2, 1}, {3, 0} },
        { {0, -1}, {1, 1}, {2, -2} },
        { {1, 0}, {0, 1}, {1, -1} }
    };
    WilsonLine correct(correct_conjugate);


    // Test complex conjugate
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ASSERT_ALMOST_EQUAL(conjugate.Element(i,j).real(), correct_conjugate[i][j].real(), 1e-10);
            ASSERT_ALMOST_EQUAL(conjugate.Element(i,j).imag(), correct_conjugate[i][j].imag(), 1e-10);
        }
    }   
}

TEST(WilsonLineTranspose)
{
    std::vector<std::vector<std::complex<double>>> matrix =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };

    WilsonLine w(matrix);
    WilsonLine transpose = w.Transpose();
    std::vector<std::vector<std::complex<double>>> correct_transpose = {
        { {1, 1}, {0, 1}, {1, 0} },
        { {2, -1}, {1, -1}, {0, -1} },
        { {3, 0}, {2, 2}, {1, 1} }
    };

    // Test transpose
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ASSERT_ALMOST_EQUAL(transpose.Element(i,j).real(), correct_transpose[i][j].real(), 1e-10);
            ASSERT_ALMOST_EQUAL(transpose.Element(i,j).imag(), correct_transpose[i][j].imag(), 1e-10);
        }
    }   
}

TEST(WilsonLineElements)
{
    std::vector<std::vector<std::complex<double>>> matrix =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };
    WilsonLine w(matrix);
    WilsonLine w2(matrix);

    w2.Set(1,2,0);
    double b = w2(1,2).real();
    b=5;
    w.Set(1,2,b);

    ASSERT_ALMOST_EQUAL(w2(1,2).real(), 0, 1e-10);

    ASSERT_ALMOST_EQUAL(w(1,2).real(), 5, 1e-10);
}

TEST(WilsonLineScalarMultiplication)
{
    std::vector<std::vector<std::complex<double>>> matrix =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };
    std::complex<double> scalar = {2, -1};

    WilsonLine w(matrix);
    WilsonLine product = w * scalar;
    std::vector<std::vector<std::complex<double>>> correct_product = {
        { {3, 1}, {3, -4}, {6, -3} },
        { {1, 2}, {1, -3}, {6, 2} },
        { {2, -1}, {-1, -2}, {3, 1} }
    };

    // Test scalar multiplication
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ASSERT_ALMOST_EQUAL(product.Element(i,j).real(), correct_product[i][j].real(), 1e-10);
            ASSERT_ALMOST_EQUAL(product.Element(i,j).imag(), correct_product[i][j].imag(), 1e-10);
        }
    }   
}


TEST(WilsonLineBounds)
{
    std::vector<std::vector<std::complex<double>>> matrix =  { {{1, 1}, {2, -1}, {3, 0}}, {{0, 1}, {1, -1}, {2, 2}}, {{1, 0}, {0, -1}, {1, 1}} };
    WilsonLine w(matrix);

    bool success=true;
    try{ 
        w.Set(4,1,1);
    }
    catch (const std::out_of_range& oor) {
        success=false;
        std::cerr << "Out of bounds error test: " << oor.what() << std::endl;
    }
    ASSERT_FALSE(success);
}

TEST_MAIN()