/*
 * GradientPSF.hpp
 *
 *  Created on: 26 mag 2017
 *      Author: Enzo
 */

#pragma once
#include <string>
#include "PointPSF.hpp"
#include "Gradient.hpp"
#include "Fourier.hpp"
#include "Maintainance.hpp"
#include <vector>

using namespace std;

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientPSF : public Gradient<T>
		{
		public:
			GradientPSF(T pixelWidth, ushort_t nDataPoints);
			virtual ~GradientPSF() { }
			virtual void updateData();
			virtual std::string getType() { return "Gradient.PSF"; }
		protected:
			T pixelWidth;
			FourierAnalysis::Fourier<T> fftR;
		};

		template <typename T>
		GradientPSF<T>::GradientPSF(T pixelWidth, ushort_t nDataPoints)
		{
			this->currentPoint = new Point::PointPSF<T>();
			this->pixelWidth = pixelWidth;
			this->currentPoint->get("sigma") = pixelWidth / GAUSS_FWHM;

			for(ushort_t x = 0; x < nDataPoints; x++)
				this->currentGradient.push_back(new Point::PointPSF<T>());

			T sigma = pixelWidth / GAUSS_FWHM;
			sigma *= sigma;
			boost::shared_ptr<T[]> R(new T[nDataPoints]);
			T flux = 0.;
			for(ushort_t x = 0; x < nDataPoints; x++)
			{
				int x2 = std::min(x * x, (nDataPoints - x) * (nDataPoints - x)); // conservazione delle posizioni
				R[x] = exp(-0.5 * x2 / sigma);
				flux += R[x];
			}
			for(ushort_t x = 0; x < nDataPoints; x++)
				R[x] /= flux; // conservazione del flusso
			this->fftR(R, nDataPoints);
		}

		template <typename T>
		void GradientPSF<T>::updateData()
		{
			boost::shared_ptr<T[]> a(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> b(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> c(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> beta(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> zero(new T[this->nDataPoints]);
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				T A = 1 + this->currentPoint->get("b") * (x - this->currentPoint->get("c")) *
						(x - this->currentPoint->get("c"));
				T B = exp(-1 * this->currentPoint->get("beta") * log(A));
				T C = B / A;
				a[x] = B;
				b[x] = -1 * this->currentPoint->get("amp") * this->currentPoint->get("beta") *
						(x - this->currentPoint->get("c")) * (x - this->currentPoint->get("c")) * C;
				c[x] = 2 * this->currentPoint->get("amp") * this->currentPoint->get("b") *
						this->currentPoint->get("beta") * (x - this->currentPoint->get("c")) * C;
				beta[x] = -1 * this->currentPoint->get("amp") * B * log(A);
				zero[x] = 1.;
			}
			this->fftR.ConvolveWith(a, this->nDataPoints, a, NULL);
			this->fftR.ConvolveWith(b, this->nDataPoints, b, NULL);
			this->fftR.ConvolveWith(c, this->nDataPoints, c, NULL);
			this->fftR.ConvolveWith(beta, this->nDataPoints, beta, NULL);
			this->fftR.ConvolveWith(zero, this->nDataPoints, zero, NULL);
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				this->currentGradient[x]["amplitude"] = a[x];
				this->currentGradient[x]["b"] = b[x];
				this->currentGradient[x]["center"] = c[x];
				this->currentGradient[x]["beta"] = beta[x];
				this->currentGradient[x]["sigma"] = 0.;
				this->currentGradient[x]["const"] = zero[x];
			}
		}
	}
}
