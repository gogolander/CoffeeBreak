/*
 * FunctionMoffat.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "PointMoffat.hpp"
#include "Function.hpp"
#include <string>
#include <vector>

using namespace std;

namespace Library
{
    namespace Function
    {
        template <typename T>
        class FunctionMoffat : public Function<T>
        {
        public:
            FunctionMoffat(ushort_t nDataPoints);
            virtual ~FunctionMoffat() { }
            virtual void updateData();
            virtual string getType() { return string("Function.Moffat"); }
        };

        template <typename T>
        FunctionMoffat<T>::FunctionMoffat(ushort_t nDataPoints)
        {
            this->nDataPoints = nDataPoints;
            this->currentData = vector<T>(nDataPoints);
        }

		template <typename T>
		void FunctionMoffat<T>::updateData()
		{
			for(ushort_t x = 0; x < this->nDataPoints; x++)
			{
				T amplitude = this->currentPoint->get("amplitude");
				T xc2 = static_cast<T>(x) - this->currentPoint->get("center");
				xc2 *= xc2;
				this->currentData[x] = amplitude * pow(1 + this->currentPoint->get("b") * xc2,
						-1 * this->currentPoint->get("beta")) +
						this->currentPoint->get("const");
			}
		}
	}
}
