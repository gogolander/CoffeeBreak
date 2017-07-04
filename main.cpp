/*
 * main.cpp
 *
 *  Created on: 09 mag 2017
 *      Author: Enzo
 */
#include <iostream>
#include <string>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>
#include <complex>
#include "Fourier.hpp"
#include "FitGauss.hpp"
#include "FitMoffat.hpp"
//#include "FitAGN.hpp"
//#include "FitPSF.hpp"
#include <vector>
//#include "Point.hpp"
#include "PointGauss.hpp"
#include "PointMoffat.hpp"
//#include "PointPSF.hpp"
//#include "PointAGN.hpp"
//#include "Function.hpp"
#include "FunctionGauss.hpp"
#include "ImageFactory.hpp"
#include <cmath>

using namespace std;
using namespace Library::Point;
using namespace Library::Function;
using namespace Library::Fit;

int main()
{
    typedef double mytype;
	PointGauss<mytype> point;
	point["amp"] = 10;
	FunctionGauss<mytype> f(64);
	point["const"] = 0.;
	typedef boost::multi_array<mytype, 2> matrix;
	matrix data(boost::extents[2048][64]);
	cout << "Crunching numbers..." << endl;
	for(matrix::index row = 0; row < 2048; row++)
	{
		point["amp"] = 10000.;
		point["center"] = 32. + 5. * sin(row/ 2048. * 6.28);
		point["sigma"] = 5. * (0.1 * row / 2048. + 1.0);
		vector<mytype> dataRow(64);
		for(matrix::index col = 0; col < 64; col++)
		{
			data[row][col] = f[point][col];
			dataRow[col] = f[point][col];
		}
		FitGauss<mytype> prova(dataRow, 64);
		PointGauss<mytype> firstGuess;
		firstGuess["amplitude"] = 9000.;
		firstGuess["center"] = 32.;
		firstGuess["sigma"] = 1.0;
//		firstGuess["beta"] = 4.5;
		firstGuess["zero"] = 0.;
//		cout << "First guess: " << (string)firstGuess << endl;
		prova.LevenbergMarquardt(firstGuess, 1E-8, 1E-8, 1E-8, 1000);
		cout << "Final point: " << (string)firstGuess << endl;
	}
	Library::ImageFactory<mytype>::writeImage("prova.fits", data, 2048, 64, true, "", true);
	cout << "Yeeaaaayyy!" << endl << "Team Alpha: We're movin' out!" << endl;
	return 0;
}
