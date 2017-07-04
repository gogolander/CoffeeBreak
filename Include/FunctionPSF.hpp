/*
 * FunctionPSF.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointPSF.hpp"
#include "Function.hpp"
#include "Fourier.hpp"
#include "Maintainance.hpp"
#include <cmath>

namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionPSF : public Function<T>
		{
		public:
			FunctionPSF(T pixelWidth, ushort_t nDataPoints);
			void updateData();
			std::string getType() { return "Function.PSF"; }
		private:
			T pixelWidth;
			FourierAnalysis::Fourier<T> fftR;
		};

		template <typename T>
		FunctionPSF<T>::FunctionPSF(T pixelWidth, ushort_t nDataPoints)
		{
			this->currentPoint = new Point::PointPSF<T>();
			this->currentData = vector<T>(nDataPoints);
			this->pixelWidth = pixelWidth;
			this->currentPoint->get("sigma") = pixelWidth / GAUSS_FWHM;
			vector<T> R(nDataPoints);
			T flux = 0.;
			for(ushort_t x = 0; x < nDataPoints; x++)
			{
				int x2 = std::min(x * x, (nDataPoints - x) * (nDataPoints - x)); // conservazione delle posizioni
				R[x] = exp(-0.5 * x2 / (this->currentPoint->get("sigma") * this->currentPoint->get("sigma")));
				flux += R[x];
			}
			for(ushort_t x = 0; x < nDataPoints; x++)
				R[x] /= flux; // conservazione del flusso
			this->fftR(R, nDataPoints);
		}

		template <typename T>
		void FunctionPSF<T>::updateData()
		{
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				T amplitude = this->currentPoint->get("amplitude");
				T xc2 = static_cast<T>(x) - this->currentPoint->get("center");
				xc2 *= xc2;
				this->currentData[x] = amplitude * pow(1 + this->currentPoint->get("b") * xc2,
						-1 * this->currentPoint->get("beta")) + this->currentPoint->get("const");
			}
			this->fftR.ConvolveWith(this->currentData, this->nDataPoints, this->currentData, NULL);
		}
	}
}
