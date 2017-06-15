/*
 * GradientMoffat.hpp
 *
 *  Created on: 26 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointMoffat.hpp"
#include "Gradient.hpp"

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientMoffat : public Gradient<T>
		{
		public:
			GradientMoffat(Index N);
			virtual ~GradientMoffat() { }
			virtual void updateData();
			virtual std::string getType() { return "Gradient.Moffat"; }
		};

		template <typename T>
		GradientMoffat<T>::GradientMoffat(Index N) : Gradient<T>(N)
		{
			this->currentPoint = new Point::PointMoffat<T>();
			this->currentData.reserve(N);
			for(Index i = 0; i < N; i++)
				this->currentData[i] = new Point::PointMoffat<T>();
		}

		template <typename T>
		void GradientMoffat<T>::updateData()
		{
			for(Index x = 0; x < this->N; x++)
			{
				T A = 1 + this->currentPoint->operator[]("b") * (x - this->currentPoint->operator[]("c")) *
						(x - this->currentPoint->operator[]("c"));
				T B = exp(-1 * this->currentPoint->operator[]("beta") * log(A));
				T C = B / A;
				this->currentData[x]->operator[]("amplitude") = B;
				this->currentData[x]->operator[]("b") = this->currentPoint->operator[]("a") *
					this->currentPoint->operator[]("beta") * (x - this->currentPoint->operator[]("c")) *
					(x + this->currentPoint->operator[]("c")) * C;
				this->currentData[x]->operator[]("center") = 2 * this->currentPoint->operator[]("a") *
						this->currentPoint->operator[]("b") * this->currentPoint->operator[]("beta") *
						(x - this->currentPoint->operator[]("c")) * C;
				this->currentData[x]->operator[]("beta") = -1 * this->currentPoint->operator[]("a") * B * log(A);
				this->currentData[x]->operator[]("const") = 1.;
			}
		}
	}
}
