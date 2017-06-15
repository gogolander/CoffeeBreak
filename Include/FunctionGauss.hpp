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
namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionGauss : public Function<T>
		{
		public:
			FunctionGauss(Index N);
			void updateData();
			Vector<T>& operator[](Point::Point<T>& other);
			std::string getType() { return "Function.Gauss"; }
		};

		template <typename T>
		FunctionGauss<T>::FunctionGauss(Index N) : Function<T>(N)
		{
		}

		template <typename T>
		Vector<T>& FunctionGauss<T>::operator[](Point::Point<T>& other)
		{
			if(!this->currentPoint)
			{
				this->currentPoint = new Point::PointGauss<T>();
				other.copyTo(this->currentPoint);
				updateData();
			}
			else if(!this->currentPoint->equals(other))
			{
				other.copyTo(this->currentPoint);
				updateData();
			}
			return this->currentData;
		}

		template <typename T>
		void FunctionGauss<T>::updateData()
		{
			if(this->currentPoint != 0)
			{
				this->currentData.clear();
				for(Index x = 0; x < this->N; x++)
				{
					T amplitude = this->currentPoint->get("amplitude");
					T fac1 = x - this->currentPoint->get("center");
					fac1 *= fac1;
					T fac2 = this->currentPoint->get("sigma");
					fac2 *= fac2;
					this->currentData.push_back(amplitude * exp(-0.5 * fac1 / fac2) + this->currentPoint->get("const"));
				}
			}
		}
	}
}
