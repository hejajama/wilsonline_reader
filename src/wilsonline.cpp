/*
 * Diffraction at sub-nucleon scale
 * Wilson line / SU(3) matrix handling
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2015-2025
 * Standalone class, no external dependences, please!
 */

#include "wilsonline.hpp"
#include <vector>
#include <complex>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;




void WilsonLine::Set(int row, int column, std::complex<double> value)
{
    if (row >= NC)
    {
        cerr << "Invalid row index " << row << " num of rows in matrix " << NC << endl;
        return;
    }
    if (column >= NC )
    {
        cerr << "Invalid column index " << column << " num of rows in matrix " << NC << endl;
        return;
    }

    operator()(row,column) = value;

}

WilsonLine WilsonLine::operator*(WilsonLine& w)
{

    WilsonLine result;

    for (int i = 0; i < NC; ++i) {
        for (int j = 0; j < NC; ++j) {
            std::complex<double> sum(0.0, 0.0);
            for (int k = 0; k < NC; ++k) {
                sum += operator()(i, k) * w(k, j);
            }
            result(i, j)=sum;
        }
    }


    return result;
}

// Multiplies this by w^\dagger, returns the product
// Fast in BLAS
WilsonLine WilsonLine::MultiplyByHermitianConjugate(const WilsonLine& w)
{
    WilsonLine result;
    for (int i = 0; i < NC; ++i) {
        for (int j = 0; j < NC; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < NC; ++k) {
                result(i, j) += (*this)(i, k) * std::conj(w(j, k));
            }
        }
    }
    return result;
    
}

std::complex<double> WilsonLine::Trace() const {
    std::complex<double> trace = 0;
    for (int i = 0; i < NC; ++i) {
        trace += (*this)(i, i);
    }
    return trace;
}

WilsonLine WilsonLine::operator*(std::complex<double> t)
{
    std::cerr << "Multiplication by a number not implemented" << endl;
    exit(1);
}


WilsonLine WilsonLine::operator+(WilsonLine& w)
{
    WilsonLine result;
    for (int i = 0; i < NC; ++i) {
        for (int j = 0; j < NC; ++j) {
            result(i, j) = (*this)(i,j)+w(i,j);
        }
    }

    return result;
}




WilsonLine WilsonLine::ComplexConjugate()
{

    std::cerr << "WilsonLine:ComplexConjugate not yet implemented for GSL matrix" << endl;
    exit(1);

}

WilsonLine WilsonLine::Transpose()
{

    std::cerr << "WilsonLine:Transpose not yet implemented for GSL matrix" << endl;
    exit(1);

}

WilsonLine WilsonLine::HermitianConjugate()
{
    
    WilsonLine res = ComplexConjugate().Transpose();
    return res;
}


int WilsonLine::Size()
{
    return NC;
}

std::complex<double> WilsonLine::Element(int row , int col) const
{
    return data[row][col];

}
std::ostream& operator<<(std::ostream& os, WilsonLine& wl)
{
    for (int r=0; r<wl.Size(); r++)
    {
        for (int i=0; i<wl.Size(); i++)
            os << wl.Element(r,i) ;
        os << endl;
    }
    return os;
}

/* Initialize as SU(3) generator
 * a=1...8
 */
void WilsonLine::InitializeAsGenerator(int a)
{
    if (a < 1 or a >8)
    {
        std::cerr << "SU(3) generator id must be within 1...8, asked " << a << std::endl;
        return;
    }
    
    /*
    std::complex<double> imag(0,1.0);
    std::vector< std::complex<double> > row;
    switch(a)
    {
        case 1:
            row.clear();
            row.push_back(0); row.push_back(0.5); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 2:
            row.clear();
            row.push_back(0); row.push_back(-0.5*imag); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5*imag); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 3:
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(-0.5); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 4:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0.5);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 5:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(-0.5*imag);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5*imag); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 6:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0.5);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0.5); row.push_back(0);
            data.push_back(row);
            break;
        case 7:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(-0.5*imag);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0.5*imag); row.push_back(0);
            data.push_back(row);
            break;
        case 8:
            row.clear();
            row.push_back(1.0/(2.0*std::sqrt(3))); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(1.0/(2.0*std::sqrt(3))); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(-2.0/(2.0*std::sqrt(3)));
            data.push_back(row);
            break;
        default:
            std::cerr << "Incorrect a!" << std::endl;
    };
    */
    cerr << "Generator not yet implemented! " << endl;
    exit(1);
}

void WilsonLine::InitializeAsIdentity()
{
    
    for (int i=0; i<NC; i++)
    {
        for (int j=0; j<NC; j++)
        {
            if (i==j)
                data[i][j]=1;
            else
                data[i][j]=0;
        }
    }
}



/*
 * Calculate exp of matrix using GSL
 * Way too much memory allocation, so one should optimize if doing some
 * heavy numerics
 */
WilsonLine WilsonLine::Exp()
{
    cerr << "Matrix exponent not impelmented" << endl;
    exit(1);
    
    return WilsonLine();
}



