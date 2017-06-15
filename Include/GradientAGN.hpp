/*
 * GradientAGN.hpp
 *
 *  Created on: 26 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointAGN.hpp"
#include "Gradient.hpp"
#include "Fourier.hpp"
#include <cmath>

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientAGN : public Gradient<T>
		{
		public:
			GradientAGN(Point::PointMoffat<T>& pointS, Index N);
			virtual ~GradientAGN() {}
			virtual void updateData();
			virtual std::string getType() { return "Gradient.AGN"; }
		protected:
			FourierAnalysis::Fourier<T> fftS;
		};

		template <typename T>
		GradientAGN<T>::GradientAGN(Point::PointMoffat<T>& pointS, Index N) : Gradient<T>(N)
		{
			this->currentPoint = new Point::PointAGN<T>();
			this->currentData(new Point::PointAGN<T>[N]);
			boost::shared_ptr<T[]> S(new T[N]);
			T flux = 0.;
			for(Index x = 0; x < N; x++)
			{
				int x2 = std::min(x * x, (N - x) * (N - x)); // conservazione delle posizioni
				S[x] = pow(1 + this->pointS["b"] * x2, -1 * this->pointS["beta"]);
				flux += S[x];
			}
			for(Index x = 0; x < N; x++)
				S[x] /= flux; // conservazione del flusso
			this->fftS(S, N);
		}

		template <typename T>
		void GradientAGN<T>::updateData()
		{
			boost::shared_ptr<T[]> aR(new T[this->N]);
			boost::shared_ptr<T[]> cR(new T[this->N]);
			boost::shared_ptr<T[]> sigmaR(new T[this->N]);
			boost::shared_ptr<T[]> aH(new T[this->N]);
			boost::shared_ptr<T[]> bH(new T[this->N]);
			boost::shared_ptr<T[]> cH(new T[this->N]);
			boost::shared_ptr<T[]> betaH(new T[this->N]);
			boost::shared_ptr<T[]> zero(new T[this->N]);
			for(Index x = 0; x < this->N; x++)
			{
				T a1 = (x - this->currentPoint["bh center"]) / this->currentPoint["bh sigma"];
				T ex = exp(-0.5 * a1 * a1);
				T fac = 2. * this->currentPoint["bh amp"] * ex * a1;
				aR[x] = ex;
				cR[x] = fac / this->currentPoint["bh sigma"];
				sigmaR = fac * a1 / this->currentPoint["bh sigma"];
				T AH = 1 + this->currentPoint["host b"] * x * x -
						2 * this->currentPoint["host b"] * x * this->currentPoint["host b"] +
						this->currentPoint["host b"] * this->currentPoint["host c"] * this->currentPoint["host c"];
				T BH = exp(-1 * this->currentPoint["host beta"] * log(AH));
				T CH = BH / AH;
				aH[x] = BH;
				bH[x] = this->currentPoint["host a"] * this->currentPoint["host beta"] *
						(2 * x * this->currentPoint["host c"] - x * x - this->currentPoint["host c"] * this->currentPoint["host c"]) * CH;
				cH[x] = 2 * this->currentPoint["host a"] * this->currentPoint["host b"] *
						this->currentPoint["host beta"] * (x - this->currentPoint["host c"]) * CH;
				betaH[x] = -1 * this->currentPoint["host a"] * BH * log(AH);
				zero[x] = 1.;
			}
			this->fftS.ConvolveWith(aR, this->N, aR, NULL);
			this->fftS.ConvolveWith(cR, this->N, cR, NULL);
			this->fftS.ConvolveWith(sigmaR, this->N, sigmaR, NULL);
			this->fftS.ConvolveWith(aH, this->N, aH, NULL);
			this->fftS.ConvolveWith(bH, this->N, bH, NULL);
			this->fftS.ConvolveWith(cH, this->N, cH, NULL);
			this->fftS.ConvolveWith(betaH, this->N, betaH, NULL);
			this->fftS.ConvolveWith(zero, this->N, zero, NULL);
			for(Index x = 0; x < this->N; x++)
			{
				this->currentData[x]["bh amplitude"] = aR[x];
				this->currentData[x]["bh center"] = cR[x];
				this->currentData[x]["bh sigma"] = sigmaR[x];
				this->currentData[x]["host amplitude"] = aH[x];
				this->currentData[x]["host b"] = bH[x];
				this->currentData[x]["host center"] = cH[x];
				this->currentData[x]["host beta"] = betaH[x];
				this->currentData[x]["const"] = zero[x];
			}
		}
	}
}
