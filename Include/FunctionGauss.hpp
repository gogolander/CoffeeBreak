/*
 * FunctionGauss.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include <cmath>
#include "PointGauss.hpp"
#include "Function.hpp"
#include <string>
#include <vector>

using namespace std;

namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionGauss : public Function<T>
		{
		public:
			FunctionGauss(ushort_t N);
			virtual ~FunctionGauss() { }
			virtual void updateData();
			virtual string getType() { return "Function.Gauss"; }
		};

		template <typename T>
		FunctionGauss<T>::FunctionGauss(ushort_t nDataPoints)
		{
			this->nDataPoints = nDataPoints;
			this->currentData = vector<T>(nDataPoints);
		}

		template <typename T>
		void FunctionGauss<T>::updateData()
		{
			this->currentData.clear();
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				T amplitude = this->currentPoint->get("amplitude");
				T fac1 = x - this->currentPoint->get("center");
				fac1 *= fac1;
				T fac2 = this->currentPoint->get("sigma");
				fac2 *= fac2;
				this->currentData[x] = amplitude * exp(-0.5 * fac1 / fac2) + this->currentPoint->get("const");
			}
		}
	}
}
