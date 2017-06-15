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

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientGauss : public Gradient<T>
		{
		public:
			GradientGauss(Index N);
			virtual ~GradientGauss() { }
			virtual std::vector<Point::Point<T>* >& operator[](Point::Point<T>& other);
			virtual void updateData();
			virtual std::string getType() { return "Gradient.Gauss"; }
		};

		template <typename T>
		GradientGauss<T>::GradientGauss(Index N) : Gradient<T>(N)
		{
			this->currentPoint = new Point::PointGauss<T>();
			for(Index i = 0; i < N; i++)
				this->currentData[i] = new Point::PointGauss<T>();
		}

		template <typename T>
		void GradientGauss<T>::updateData()
		{
			for(Index x = 0; x < this->N; x++)
			{
				T a1 = (x - this->currentPoint->get("center")) / this->currentPoint->get("sigma");
				T ex = exp(-0.5 * a1 * a1);
				T fac = 2. * this->currentPoint->get("amp") * ex * a1;
				this->currentData[x]->get("amplitude") = ex;
				this->currentData[x]->get("center") = fac / this->currentPoint->get("sigma");
				this->currentData[x]->get("sigma") = fac * a1 / this->currentPoint->get("sigma");
				this->currentData[x]->get("const") = 1.;
			}
		}
	}
}
