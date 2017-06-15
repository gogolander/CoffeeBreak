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
#include "FitAGN.hpp"
#include "FitPSF.hpp"
#include <vector>
#include "Point.hpp"
#include "PointGauss.hpp"
#include "PointMoffat.hpp"
#include "PointPSF.hpp"
#include "PointAGN.hpp"
#include "Function.hpp"
#include "FunctionGauss.hpp"
#include "ImageFactory.hpp"

using namespace std;
using namespace Library::Point;
using namespace Library::Function;

int main()
{
	std::size_t cols = 0;
	std::size_t rows = 0;
	std::string filename = "C:/Users/Enzo/CoffeeBreak/Resources/15030208/R/star11.fits";
	cout << rows << " rows; " << cols << " cols" << endl;
	cout << "Yeeaaaayyy!" << endl << "Team Alpha: We're movin' out!" << endl;
	PointGauss<float> point;
	point["amp"] = 10;
	PointGauss<float> clone;
	cout << (string)point << endl;
	cout << (string)clone << endl;
	clone = point;
	cout << (string)clone << endl;
	cout << (string)(clone + point) << endl;
	FunctionGauss<float> f(64);
//	point["const"] = 0.;
//	typedef boost::multi_array<float, 2> matrix;
//	matrix data(boost::extents[2048][64]);
//	for(matrix::index row = 0; row < 2048; row++)
//	{
//		point["amp"] = row;
//		point["center"] = 32;
//		point["sigma"] = (row + 1.) / 2048. * 5.;
//		for(matrix::index col = 0; col < 64; col++)
//			data[row][col] = f[point][col];
//	}
//	Library::ImageFactory<float>::writeImage("prova.fits", data, 2048, 64, true, "", true);
	return 0;
}
