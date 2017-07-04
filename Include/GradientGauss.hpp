/*
 * GradientGauss.hpp
 *
 *  Created on: 26 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointGauss.hpp"
#include "Point.hpp"
#include "Gradient.hpp"
#include <vector>
#include <string>

using namespace std;

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientGauss : public Gradient<T>
		{
		public:
			GradientGauss(ushort_t nDataPoints);
			virtual ~GradientGauss() { }
			virtual void updateData();
			virtual string getType() { return string("Gradient.Gauss"); }
		};

		template <typename T>
		GradientGauss<T>::GradientGauss(ushort_t nDataPoints)
		{
			this->nDataPoints = nDataPoints;
			for(ushort_t x = 0; x < nDataPoints; x++)
				this->currentGradient.push_back(new Point::PointGauss<T>());
		}

		template <typename T>
		void GradientGauss<T>::updateData()
		{
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				T a1 = (static_cast<T>(x) - this->currentPoint->get("center")) / this->currentPoint->get("sigma");
				T ex = exp(-0.5 * a1 * a1);
				T fac = 2. * this->currentPoint->get("amp") * ex * a1;
				this->currentGradient[x]->get("amplitude") = ex;
				this->currentGradient[x]->get("center") = fac / this->currentPoint->get("sigma");
				this->currentGradient[x]->get("sigma") = fac * a1 / this->currentPoint->get("sigma");
				this->currentGradient[x]->get("const") = 1.;
			}
		}
	}
}
