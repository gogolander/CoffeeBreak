/*
 * GradientMoffat.hpp
 *
 *  Created on: 26 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointMoffat.hpp"
#include "Gradient.hpp"
#include "Maintainance.hpp"

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientMoffat : public Gradient<T>
		{
		public:
			GradientMoffat(ushort_t nDataPoints);
			virtual ~GradientMoffat() { }
			virtual void updateData();
			virtual std::string getType() { return "Gradient.Moffat"; }
		};

		template <typename T>
		GradientMoffat<T>::GradientMoffat(ushort_t nDataPoints)
		{
			this->nDataPoints = nDataPoints;
			for(ushort_t x = 0; x < nDataPoints; x++)
				this->currentGradient.push_back( new Point::PointMoffat<T>() );
		}

		template <typename T>
		void GradientMoffat<T>::updateData()
		{
			for(Index x = 0; x < this->N; x++)
			{
				T A = 1 + this->currentPoint->get("b") * (static_cast<T>(x) - this->currentPoint->get("c")) *
						(static_cast<T>(x) - this->currentPoint->get("c"));
				T B = exp(-1 * this->currentPoint->get("beta") * log(A));
				T C = B / A;
				this->currentData[x]->get("amplitude") = B;
				this->currentData[x]->get("b") = this->currentPoint->get("a") *
					this->currentPoint->get("beta") * (static_cast<T>(x) - this->currentPoint->get("c")) *
					(static_cast<T>(x) + this->currentPoint->get("c")) * C;
				this->currentData[x]->get("center") = 2 * this->currentPoint->get("a") *
						this->currentPoint->get("b") * this->currentPoint->get("beta") *
						(static_cast<T>(x) - this->currentPoint->get("c")) * C;
				this->currentData[x]->get("beta") = -1 * this->currentPoint->get("a") * B * log(A);
				this->currentData[x]->get("const") = 1.;
			}
		}
	}
}
