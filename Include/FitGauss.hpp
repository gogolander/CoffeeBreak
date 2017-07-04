/*
 * FitGauss.hpp
 *
 *  Created on: 27 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Fit.h"
#include "Fit.cpp"
#include "PointGauss.hpp"
#include "FunctionGauss.hpp"
#include "GradientGauss.hpp"
#include "Maintainance.hpp"
#include "Data.hpp"

#include <iostream>
#include <boost/multi_array.hpp>
#include <vector>

namespace Library
{
	namespace Fit
	{
		template <typename T>
		class FitGauss : public Fit<T>
		{
		public:
			FitGauss(vector<T>& data, ushort_t nDataPoints);
			void HongweiGuo(Point::PointGauss<T>& Point, int iterations);
			virtual ~FitGauss() { }
		protected:
			virtual Function::Function<T>& function() { return myFunction; }
			virtual Function::Gradient<T>& gradient() {return myGradient; }
			virtual Point::Point<T>& newPoint() { return myNewPoint; }
			virtual Point::Point<T>& deltaPoint() { return myDeltaPoint; }
			virtual Point::Point<T>& h_sd() { return my_h_sd; }
			virtual Point::Point<T>& h_gn() { return my_h_gn; }
			virtual Point::Point<T>& h_dl() { return my_h_dl; }

			Function::FunctionGauss<T> myFunction;
			Function::GradientGauss<T> myGradient;
			Point::PointGauss<T> myNewPoint;
			Point::PointGauss<T> myDeltaPoint;
			Point::PointGauss<T> my_h_sd;
			Point::PointGauss<T> my_h_gn;
			Point::PointGauss<T> my_h_dl;
		};

		template<typename T>
		FitGauss<T>::FitGauss(vector<T>& data, ushort_t nDataPoints) : myFunction(nDataPoints),
			myGradient(nDataPoints)
		{
			this->data = Data::Data<T>(data, nDataPoints);
			this->nDataPoints = nDataPoints;
			this->data.initData();
		}

		template <typename T>
		void FitGauss<T>::HongweiGuo(Point::PointGauss<T>& Point, int iterations)
		{
			T** A = new T*[3];
			for(int i = 0; i < 3; i++)
				A[i] = new T[3];
			T b[3];
			T solution[3];
			A[0][0] = A[0][1] = A[0][2] = A[1][2] = A[2][2] = b[0] = b[1] = b[2] = solution[0] = solution[1] = solution[2] = 0; // these are the only critical elements to set to zero thanks to symmetry
			for(int x = 0; x < this->nDataPoints; x++)
			{
				A[0][0] += this->data[x];
				A[0][1] += x * A[0][0];
				A[0][2] += x * A[0][1];
				A[1][2] += x * A[0][2];
				A[2][2] += x * A[1][2];
				b[0] += this->data[x] * this->data[x] * log(this->data[x]);
				b[1] += x * b[0];
				b[2] += x * b[1];
			}
			A[1][0] = A[0][1];
			A[2][0] = A[1][1] = A[0][2];
			A[2][1] = A[1][2];
			this->Inverse(A, 3);
			for(int i = 0; i < 3; i++)
				for(int j = 0; j < 3; j++)
					solution[i] += A[j][i] * b[j];
			for(int k = 1; k < iterations; k++)
			{
				A[0][0] = A[0][1] = A[0][2] = A[1][2] = A[2][2] = b[0] = b[1] = b[2] = 0; // these are the only critical elements to set to zero thanks to symmetry
				for(int x = 0; x < this->nDataPoints; x++)
				{
					T y = exp(solution[0] + x * solution[1] + x * x * solution[2]);
					A[0][0] += y;
					A[0][1] += x * A[0][0];
					A[0][2] += x * A[0][1];
					A[1][2] += x * A[0][2];
					A[2][2] += x * A[1][2];
					b[0] += y * y * log(y);
					b[1] += x * b[0];
					b[2] += x * b[1];
				}
				A[1][0] = A[0][1];
				A[2][0] = A[1][1] = A[0][2];
				A[2][1] = A[1][2];
				this->Inverse(A, 3);
				for(int i = 0; i < 3; i++)
				{
					solution[i] = 0;
					for(int j = 0; j < 3; j++)
						solution[i] += A[j][i] * b[j];
				}
			}
			Point[0] = exp(solution[0] - solution[1] * solution[1] / (4 * solution[2]));
			Point[1] = -1 * solution[1] / (2 * solution[2]);
			Point[2] = sqrt(-1 / (2 * solution[2]));
			delete[] A;
		}
	}
}
