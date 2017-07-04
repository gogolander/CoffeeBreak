/*
 * FitAGN.hpp
 *
 *  Created on: 27 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Fit.h"
#include "PointAGN.hpp"
#include "PointMoffat.hpp"
#include "FunctionAGN.hpp"
#include "GradientAGN.hpp"

namespace Library {
    namespace Fit {
        template<typename T>
        class FitAGN: public Fit<T> {
            public:
                FitAGN(vector<T>& data, ushort_t nDataPoints,
                        Point::PointMoffat<T>& pointS);
                virtual ~FitAGN() {
                }
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

                Function::FunctionAGN<T> myFunction;
                Function::GradientAGN<T> myGradient;
                Point::PointAGN<T> myNewPoint;
                Point::PointAGN<T> myDeltaPoint;
                Point::PointAGN<T> my_h_sd;
                Point::PointAGN<T> my_h_gn;
                Point::PointAGN<T> my_h_dl;
        };

        template<typename T>
        FitAGN<T>::FitAGN(vector<T>& data, ushort_t nDataPoints,
                Point::PointMoffat<T>& pointS) :
                myFunction(nDataPoints), myGradient(nDataPoints) {
            this->data = Data::Data<T>(data, nDataPoints);
            this->nDataPoints = nDataPoints;
            this->data.initData();
        }
    }
}
