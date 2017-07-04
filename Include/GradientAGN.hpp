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
#include <string>

using namespace std;

namespace Library
{
	namespace Function
	{
		template <typename T>
		class GradientAGN : public Gradient<T>
		{
		public:
			GradientAGN(Point::PointMoffat<T>& pointS, ushort_t nDataPoints);
			virtual ~GradientAGN() {}
			virtual void updateData();
			virtual string getType() { return string("Gradient.AGN"); }
		protected:
			FourierAnalysis::Fourier<T> fftS;
		};

		template <typename T>
		GradientAGN<T>::GradientAGN(Point::PointMoffat<T>& pointS, ushort_t nDataPoints)
		{
			this->currentPoint = new Point::PointAGN<T>();
			this->nDataPoints = nDataPoints;
			for(ushort_t x = 0; x < nDataPoints; x++)
				this->currentGradient.push_back(new Point::PointAGN<T>());
			boost::shared_ptr<T[]> S(new T[nDataPoints]);
			T flux = 0.;
			for(ushort_t x = 0; x < nDataPoints; x++)
			{
				int x2 = std::min(x * x, (nDataPoints - x) * (nDataPoints - x)); // conservazione delle posizioni
				S[x] = pow(1 + this->pointS["b"] * x2, -1 * this->pointS["beta"]);
				flux += S[x];
			}
			for(ushort_t x = 0; x < nDataPoints; x++)
				S[x] /= flux; // conservazione del flusso
			this->fftS(S, nDataPoints);
		}

		template <typename T>
		void GradientAGN<T>::updateData()
		{
			boost::shared_ptr<T[]> aR(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> cR(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> sigmaR(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> aH(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> bH(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> cH(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> betaH(new T[this->nDataPoints]);
			boost::shared_ptr<T[]> zero(new T[this->nDataPoints]);
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				T a1 = (x - this->currentPoint->get("bh center")) / this->currentPoint->get("bh sigma");
				T ex = exp(-0.5 * a1 * a1);
				T fac = 2. * this->currentPoint->get("bh amp") * ex * a1;
				aR[x] = ex;
				cR[x] = fac / this->currentPoint->get("bh sigma");
				sigmaR = fac * a1 / this->currentPoint->get("bh sigma");
				T AH = 1 + this->currentPoint->get("host b") * (x - this->currentPoint->get("host c")) *
						(x - this->currentPoint->get("host c"));
				T BH = exp(-1 * this->currentPoint->get("host beta") * log(AH));
				T CH = BH / AH;
				aH[x] = BH;
				bH[x] = -1* this->currentPoint->get("host a") * this->currentPoint->get("host beta") *
						(x - this->currentPoint->get("host c")) * (x - this->currentPoint->get("host c")) * CH;
				cH[x] = 2 * this->currentPoint->get("host a") * this->currentPoint->get("host b") *
						this->currentPoint->get("host beta") *
						(x - this->currentPoint->get("host c")) * CH;
				betaH[x] = -1 * this->currentPoint->get("host amp") * BH * log(AH);
				zero[x] = 1.;
			}
			this->fftS.ConvolveWith(aR, this->nDataPoints, aR, NULL);
			this->fftS.ConvolveWith(cR, this->nDataPoints, cR, NULL);
			this->fftS.ConvolveWith(sigmaR, this->nDataPoints, sigmaR, NULL);
			this->fftS.ConvolveWith(aH, this->nDataPoints, aH, NULL);
			this->fftS.ConvolveWith(bH, this->nDataPoints, bH, NULL);
			this->fftS.ConvolveWith(cH, this->nDataPoints, cH, NULL);
			this->fftS.ConvolveWith(betaH, this->nDataPoints, betaH, NULL);
			this->fftS.ConvolveWith(zero, this->nDataPoints, zero, NULL);
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				this->currentGradient[x]["bh amplitude"] = aR[x];
				this->currentGradient[x]["bh center"] = cR[x];
				this->currentGradient[x]["bh sigma"] = sigmaR[x];
				this->currentGradient[x]["host amplitude"] = aH[x];
				this->currentGradient[x]["host b"] = bH[x];
				this->currentGradient[x]["host center"] = cH[x];
				this->currentGradient[x]["host beta"] = betaH[x];
				this->currentGradient[x]["const"] = zero[x];
			}
		}
	}
}
