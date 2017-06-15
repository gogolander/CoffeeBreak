/*
 * FitMoffat.hpp
 *
 *  Created on: 27 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Fit.h"
#include "PointMoffat.hpp"
#include "FunctionMoffat.hpp"
#include "GradientMoffat.hpp"

namespace Library
{
	namespace Fit
	{
		template <typename T>
		class FitMoffat : public Fit<T>
		{
		public:
			typedef typename boost::shared_ptr<T**> Matrix;
			FitMoffat(Vector<T>& data, int nDataPoints);
		};

		template <typename T>
		FitMoffat<T>::FitMoffat(Vector<T>& data, int nDataPoints) : Fit<T>(data, nDataPoints)
		{
			this->function = new Function::FunctionMoffat<T>(nDataPoints);
			this->gradient = new Function::GradientMoffat<T>(nDataPoints);
		}
	}
}
