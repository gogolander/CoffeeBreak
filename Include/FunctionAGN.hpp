/*
 * FunctionAGN.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointAGN.hpp"
#include "Function.hpp"
#include "Fourier.hpp"
#include "Maintainance.hpp"
namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionAGN : public Function<T>
		{
		public:
			FunctionAGN(Point::PointMoffat<T>& pointS, Index N);
			Vector<T>& operator[](Point::Point<T>& other);
			void updateData();
			std::string getType() { return "Function.AGN"; }
		private:
			FourierAnalysis::Fourier<T> fftS;
		};

		template <typename T>
		FunctionAGN<T>::FunctionAGN(Point::PointMoffat<T>& pointS, Index N) : Function<T>(N)
		{
			this->currentPoint(new Point::PointAGN<T>());
			boost::shared_ptr<T[]> S(new T[N]);
			T flux = 0.;
			for(int x = 0; x < N; x++)
			{
				int x2 = std::min(x * x, (N - x) * (N - x)); // conservazione delle posizioni
				S[x] = pow(1 + this->pointS["b"] * x2, -1 * this->pointS["beta"]);
				flux += S[x];
			}
			for(int x = 0; x < N; x++)
				S[x] /= flux; // conservazione del flusso
			this->fftS(S, N);
		}

		template <typename T>
		Vector<T>& FunctionAGN<T>::operator[](Point::Point<T>& other)
		{
			if(!other.equals(this->currentPoint))
			{
				other.copyTo(this->currentPoint);
				updateData();
			}
			return this->currentData;
		}

		template <typename T>
		void FunctionAGN<T>::updateData()
		{
			for(Index x = 0; x < this->N; x++)
			{
				T fac1 = x - *(this->currentPoint.get())["bh center"];
				fac1 *= fac1;
				T fac2 = *(this->currentPoint.get())["bh sigma"];
				fac2 *= -2 * fac2;
				T bh = *(this->currentPoint.get())["bh amplitude"] * exp(fac1 / fac2);
				T hostXC2 = x - *(this->currentPoint.get())["host center"];
				hostXC2 *= hostXC2;
				T host = *(this->currentPoint.get())["host amplitude"] * pow(1 +
						*(this->currentPoint.get())["host b"] * hostXC2, -1 * *(this->currentPoint.get())["host beta"]);
				this->currentData[x] = bh + host + *(this->currentPoint.get())["const"];
			}
			this->fftR.ConvolveWith(this->currentData, this->N, this->currentData, NULL);
		}
	}
}
