/*
 * FitPSF.hpp
 *
 *  Created on: 27 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Fit.h"
#include "PointPSF.hpp"
#include "PointGauss.hpp"
#include "FunctionPSF.hpp"
#include "GradientPSF.hpp"

namespace Library {
    namespace Fit {
        template<typename T>
        class FitPSF: public Fit<T> {
            public:
                FitPSF(vector<T>& data, ushort_t nDataPoints,
                        Point::PointGauss<T>& pointR);
                virtual ~FitPSF() {
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

                Function::FunctionPSF<T> myFunction;
                Function::GradientPSF<T> myGradient;
                Point::PointPSF<T> myNewPoint;
                Point::PointPSF<T> myDeltaPoint;
                Point::PointPSF<T> my_h_sd;
                Point::PointPSF<T> my_h_gn;
                Point::PointPSF<T> my_h_dl;
        };

        template<typename T>
        FitPSF<T>::FitPSF(vector<T>& data, ushort_t nDataPoints,
                Point::PointGauss<T>& pointR) :
                myFunction(nDataPoints), myGradient(nDataPoints) {
            this->data = Data::Data<T>(data, nDataPoints);
            this->nDataPoints = nDataPoints;
            this->data.initData();
        }
    }
}
