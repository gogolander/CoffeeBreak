/*
 * FitAGN.hpp
 *
 *  Created on: 27 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Fit.h"
#include "PointAGN.hpp"
#include "FunctionAGN.hpp"
#include "GradientAGN.hpp"

namespace Library
{
	namespace Fit
	{
		template <typename T>
		class FitAGN : public Fit<T>
		{
		public:
			typedef typename boost::shared_ptr<T**> Matrix;
			FitAGN(Vector<T>& data, int nDataPoints, Point::PointMoffat<T>& pointS);
		};

		template <typename T>
		FitAGN<T>::FitAGN(Vector<T>& data, int nDataPoints, Point::PointMoffat<T>& pointS) :
			Fit<T>(data, nDataPoints)
		{
			this->function = new Function::FunctionAGN<T>(pointS, nDataPoints);
			this->gradient = new Function::GradientAGN<T>(pointS, nDataPoints);
		}
	}
}
