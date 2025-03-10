/*
 * Diffraction at sub-nucleon scale
 * Wilson line / SU(3) matrix handling
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 * Standalone class, no external dependences, please!
 */

#ifndef wilsonline_hpp
#define wilsonline_hpp

#include <iostream>
#include <complex>
#include <vector>



#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

typedef unsigned int uint;

using std::cout;
using std::cerr;
using std::endl;




class WilsonLine 
{
public:

    WilsonLine(const std::vector<std::vector<std::complex<double>>>& mat) {
        for (int i = 0; i < NC; ++i) {
            for (int j = 0; j < NC; ++j) {
                data[i][j] = mat[i][j];
            }
        }
    }

    WilsonLine() {
        for (int i = 0; i < NC; ++i) {
            for (int j = 0; j < NC; ++j) {
                data[i][j] = 0;
            }
        }
    }

    std::complex<double>& operator()(int row, int col) {
        return data[row][col];
    }

    const std::complex<double>& operator()(int row, int col) const {
        return data[row][col];
    }

    WilsonLine MultiplyByHermitianConjugate(const WilsonLine other) const;

    std::complex<double> Trace() const; 

    WilsonLine operator*(WilsonLine& w);
    WilsonLine operator*(std::complex<double> t);
    WilsonLine operator+(WilsonLine& w);
    

    
    // Multiplies this by w^\dagger, returns the product
    WilsonLine MultiplyByHermitianConjugate(const WilsonLine& other);
    
    
    WilsonLine ComplexConjugate();
    WilsonLine Transpose();
    WilsonLine HermitianConjugate();
    
    void Set(int row, int column, std::complex<double> value);
    
    int Size(); // Size of NxN matrix
    std::complex<double> Element(int row, int col) const;
    
    void InitializeAsGenerator(int a);  // Initialize as Color matrix t^a
    void InitializeAsIdentity();
    
    WilsonLine Exp();   // Calculate exponential
private:
    static const int NC = 3;

    std::complex<double> data[NC][NC];

    

    
};


std::ostream& operator<<(std::ostream& os, WilsonLine& wl);





#endif /* wilsonline_hpp */
