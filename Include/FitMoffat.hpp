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

namespace Library {
    namespace Fit {
        template<typename T>
        class FitMoffat: public Fit<T> {
            public:
                FitMoffat(vector<T>& data, ushort_t nDataPoints);
                virtual ~FitMoffat() { }
            protected:
                virtual Function::Function<T>& function() {
                    return myFunction;
                }
                virtual Function::Gradient<T>& gradient() {
                    return myGradient;
                }
                virtual Point::Point<T>& newPoint() {
                    return myNewPoint;
                }
                virtual Point::Point<T>& deltaPoint() {
                    return myDeltaPoint;
                }
                virtual Point::Point<T>& h_sd() {
                    return my_h_sd;
                }
                virtual Point::Point<T>& h_gn() {
                    return my_h_gn;
                }
                virtual Point::Point<T>& h_dl() {
                    return my_h_dl;
                }

                Function::FunctionMoffat<T> myFunction;
                Function::GradientMoffat<T> myGradient;
                Point::PointMoffat<T> myNewPoint;
                Point::PointMoffat<T> myDeltaPoint;
                Point::PointMoffat<T> my_h_sd;
                Point::PointMoffat<T> my_h_gn;
                Point::PointMoffat<T> my_h_dl;
        };

        template<typename T>
        FitMoffat<T>::FitMoffat(vector<T>& data, ushort_t nDataPoints) :
                myFunction(nDataPoints), myGradient(nDataPoints) {
            this->data = Data::Data<T>(data, nDataPoints);
            this->nDataPoints = nDataPoints;
            this->data.initData();
        }
    }
}
