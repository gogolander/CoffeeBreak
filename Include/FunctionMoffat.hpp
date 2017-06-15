/*
 * FunctionMoffat.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointMoffat.hpp"
#include "Function.hpp"
namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionMoffat : public Function<T>
		{
		public:
			FunctionMoffat(Index N);
			Vector<T>& operator[](Point::Point<T>& other);
			void updateData();
			std::string getType() { return "Function.Moffat"; }
		};

		template <typename T>
		FunctionMoffat<T>::FunctionMoffat(Index N) : Function<T>(N)
		{
			this->currentPoint = new Point::PointMoffat<T>();
		}

		template <typename T>
		Vector<T>& FunctionMoffat<T>::operator[](Point::Point<T>& other)
		{
			if(!other.equals(this->currentPoint))
			{
				other.copyTo(this->currentPoint);
				updateData();
			}
			return this->currentData;
		}

		template <typename T>
		void FunctionMoffat<T>::updateData()
		{
			for(Index x = 0; x < this->N; x++)
			{
				T amplitude = (*this->currentPoint)["amplitude"];
				T xc2 = x - (*this->currentPoint)["center"];
				xc2 *= xc2;
				this->currentData[x] = amplitude * pow(1 + (*this->currentPoint)["b"] * xc2,
						-1 * (*this->currentPoint)["beta"]) +
						(*this->currentPoint)["const"];
			}
		}

//		template <typename T>
//		void FunctionMoffat<T>::gradientAt(T x, boost::shared_ptr<T[]>& result)
//		{
//			T a = p[0];
//			T b = p[1];
//			T c = p[2];
//			T beta = p[3];
//			result[0] = result[1] = result[2] = result[3] = 0;
//			result[4] = 1;
//			// calcolo costanti che vengono fuori spesso nel calcolo del gradiente per risparmiare tempo di calcolo
//			T A = 1 + b * x * x - 2 * b * x * c + b * c * c;
//			T B = exp(-1 * beta * log(A));
//			T C = B / A;
//			result[0] = B;
//			result[1] = a * beta * (2 * x * c - x * x - c * c) * C;
//			result[2] = 2 * a * b * beta * (x - c) * C;
//			result[3] = -1 * a * B * log(A);
//		}
//
//		template <typename T>
//		void FunctionMoffat<T>::gradientAt(Point::PointMoffat& p, T x, boost::shared_ptr<T[]>& result)
//		{
//			if(p.isValid())
//			{
//				T a = p[0];
//				T b = p[1];
//				T c = p[2];
//				T beta = p[3];
//				result[0] = result[1] = result[2] = result[3] = 0;
//				result[4] = 1;
//				// calcolo costanti che vengono fuori spesso nel calcolo del gradiente per risparmiare tempo di calcolo
//				T A = 1 + b * x * x - 2 * b * x * c + b * c * c;
//				T B = exp(-1 * beta * log(A));
//				T C = B / A;
//				result[0] = B;
//				result[1] = a * beta * (2 * x * c - x * x - c * c) * C;
//				result[2] = 2 * a * b * beta * (x - c) * C;
//				result[3] = -1 * a * B * log(A);
//			}
//			else result[0] = result[1] = result[2] = result[3] = result[4] = 1E6;
//		}
	}
}
