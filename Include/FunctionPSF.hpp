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

namespace Library
{
	namespace Function
	{
		template <typename T>
		class FunctionPSF : public Function<T>
		{
		public:
			FunctionPSF(T pixelWidth, Index N);
			Vector<T>& operator[](Point::Point<T>& other);
			void updateData();
			std::string getType() { return "Function.PSF"; }
		private:
			T pixelWidth;
			FourierAnalysis::Fourier<T> fftR;
		};

		template <typename T>
		FunctionPSF<T>::FunctionPSF(T pixelWidth, Index N) : Function<T>(N)
		{
			this->currentPoint = new Point::PointPSF<T>();
			this->pixelWidth = pixelWidth;
			this->currentPoint["sigma"] = pixelWidth / GAUSS_FWHM;
			boost::shared_ptr<T[]> R(new T[N]);
			T flux = 0.;
			for(Index x = 0; x < N; x++)
			{
				int x2 = std::min(x * x, (N - x) * (N - x)); // conservazione delle posizioni
				R[x] = exp(-0.5 * x2 / (this->currentPoint["sigma"] * this->currentPoint["sigma"]));
				flux += R[x];
			}
			for(Index x = 0; x < N; x++)
				R[x] /= flux; // conservazione del flusso
			this->fftR(R, N);
		}

		template <typename T>
		Vector<T>& FunctionPSF<T>::operator[](Point::Point<T>& other)
		{
			if(!other.equals(this->currentPoint))
			{
				other.copyTo(this->currentPoint);
				updateData();
			}
			return this->currentData;
		}

		template <typename T>
		void FunctionPSF<T>::updateData()
		{
			for(Index x = 0; x < this->N; x++)
			{
				T amplitude = this->currentPoint["amplitude"];
				T xc2 = x - this->currentPoint["center"];
				xc2 *= xc2;
				this->currentData[x] = amplitude * pow(1 + this->currentPoint["b"] * xc2,
						-1 * this->currentPoint["beta"]) + this->currentPoint["const"];
			}
			this->fftR.ConvolveWith(this->currentData, this->N, this->currentData, NULL);
		}
	}
}
