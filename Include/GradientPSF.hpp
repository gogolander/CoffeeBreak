/*
 * GradientPSF.hpp
 *
 *  Created on: 26 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointPSF.hpp"
#include "Gradient.hpp"
#include "Fourier.hpp"
namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientPSF : public Gradient<T>
		{
		public:
			GradientPSF(T pixelWidth, Index N);
			virtual ~GradientPSF() { }
			virtual void updateData();
			virtual std::string getType() { return "Gradient.PSF"; }
		protected:
			T pixelWidth;
			FourierAnalysis::Fourier<T> fftR;
		};

		template <typename T>
		GradientPSF<T>::GradientPSF(T pixelWidth, Index N) : Gradient<T>(N)
		{
			this->currentPoint = new Point::PointPSF<T>();
			this->currentData(new Point::PointPSF<T>[N]);
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
		void GradientPSF<T>::updateData()
		{
			boost::shared_ptr<T[]> a(new T[this->N]);
			boost::shared_ptr<T[]> b(new T[this->N]);
			boost::shared_ptr<T[]> c(new T[this->N]);
			boost::shared_ptr<T[]> beta(new T[this->N]);
			boost::shared_ptr<T[]> zero(new T[this->N]);
			for(Index x = 0; x < this->N; x++)
			{
				T A = 1 + this->currentPoint["b"] * x * x -
						2 * this->currentPoint["b"] * x * this->currentPoint["b"] +
						this->currentPoint["b"] * this->currentPoint["c"] * this->currentPoint["c"];
				T B = exp(-1 * this->currentPoint["beta"] * log(A));
				T C = B / A;
				a[x] = B;
				b[x] = this->currentPoint["a"] * this->currentPoint["beta"] *
						(2 * x * this->currentPoint["c"] - x * x - this->currentPoint["c"] * this->currentPoint["c"]) * C;
				c[x] = 2 * this->currentPoint["a"] * this->currentPoint["b"] *
						this->currentPoint["beta"] * (x - this->currentPoint["c"]) * C;
				beta[x] = -1 * this->currentPoint["a"] * B * log(A);
				zero[x] = 1.;
			}
			this->fftR.ConvolveWith(a, this->N, a, NULL);
			this->fftR.ConvolveWith(b, this->N, b, NULL);
			this->fftR.ConvolveWith(c, this->N, c, NULL);
			this->fftR.ConvolveWith(beta, this->N, beta, NULL);
			this->fftR.ConvolveWith(zero, this->N, zero, NULL);
			for(Index x = 0; x < this->N; x++)
			{
				this->currentData[x]["amplitude"] = a[x];
				this->currentData[x]["b"] = b[x];
				this->currentData[x]["center"] = c[x];
				this->currentData[x]["beta"] = beta[x];
				this->currentData[x]["sigma"] = 0.;
				this->currentData[x]["const"] = zero[x];
			}
		}
	}
}
