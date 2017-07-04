/*
 * FunctionAGN.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointAGN.hpp"
#include "PointMoffat.hpp"
#include "Function.hpp"
#include "Fourier.hpp"
#include "Maintainance.hpp"
#include <cmath>

namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionAGN : public Function<T>
		{
		public:
			FunctionAGN(Point::PointMoffat<T>& pointS, ushort_t nDataPoints);
			void updateData();
			std::string getType() { return string("Function.AGN"); }
		private:
			FourierAnalysis::Fourier<T> fftS;
		};

		template <typename T>
		FunctionAGN<T>::FunctionAGN(Point::PointMoffat<T>& pointS, ushort_t nDataPoints)
		{
			vector<T> S(nDataPoints);
			this->currentData = vector<T>(nDataPoints);
			T flux = 0.;
			for(ushort_t x = 0; x < nDataPoints; x++)
			{
				ushort_t x2 = std::min(x * x, (nDataPoints - x) * (nDataPoints - x)); // conservazione delle posizioni
				S[x] = pow(1 + this->pointS["b"] * x2, -1 * this->pointS["beta"]);
				flux += S[x];
			}
			for(int x = 0; x < nDataPoints; x++)
				S[x] /= flux; // conservazione del flusso
			this->fftS(S, nDataPoints);
		}

		template <typename T>
		void FunctionAGN<T>::updateData()
		{
			for(Index x = 0; x < this->nDataPoints; x++)
			{
				T fac1 = static_cast<T>(x) - this->currentPoint.get("bh center");
				fac1 *= fac1;
				T fac2 = this->currentPoint.get("bh sigma");
				fac2 *= -2 * fac2;
				T bh = this->currentPoint.get("bh amplitude") * exp(fac1 / fac2);
				T hostXC2 = static_cast<T>(x) - this->currentPoint.get("host center");
				hostXC2 *= hostXC2;
				T host = this->currentPoint.get("host amplitude") * pow(1 +
						this->currentPoint.get("host b") * hostXC2, -1 * this->currentPoint.get("host beta"));
				this->currentData[x] = bh + host + this->currentPoint.get("const");
			}
			this->fftR.ConvolveWith(this->currentData, this->nDataPoints, this->currentData, NULL);
		}
	}
}
