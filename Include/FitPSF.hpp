/*
 * FitPSF.hpp
 *
 *  Created on: 27 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Fit.h"
#include "PointPSF.hpp"
#include "FunctionPSF.hpp"
#include "GradientPSF.hpp"

namespace Library
{
	namespace Fit
	{
		template <typename T>
		class FitPSF : public Fit<T>
		{
		public:
			typedef typename boost::shared_ptr<T**> Matrix;
			FitPSF(Vector<T>& data, int nDataPoints, Point::PointGauss<T>& pointR);
		};

		template <typename T>
		FitPSF<T>::FitPSF(Vector<T>& data, int nDataPoints, Point::PointGauss<T>& pointR) : Fit<T>(data, nDataPoints)
		{
			this->function = new Function::FunctionPSF<T>(pointR, nDataPoints);
			this->gradient = new Function::GradientPSF<T>(pointR, nDataPoints);
		}
	}
}
