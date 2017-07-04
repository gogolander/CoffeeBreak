#pragma once
#include <string>
#include "Point.hpp"
#include "Function.hpp"
#include "Gradient.hpp"
#include "Data.hpp"
#include "Maintainance.hpp"
#include <vector>

using namespace std;

namespace Library
{
	namespace Fit
	{
		//<summay>Class that contains the tools for fitting experimental data with a function defined by the user.</summary>
		template <typename U>
		class Fit
		{
		public:
			Fit();
			virtual ~Fit() { }
			int LevenbergMarquardt(Point::Point<U>& point, U gradientTOL, U xTOL, U chi2TOL, int maxIterations);
			int DogLeg(Point::Point<U>& point, U Delta0, U gradientTOL, U xTOL, U chi2TOL, int maxIterations);
			Data::Data<U>& getData();
			void setData(vector<U>& data);
			vector<U>& getModel(Point::Point<U>& point);
			vector<U>& getResiduals(Point::Point<U>& point);
			U abs(vector<U>& vec);

		protected:
//			void SolveLinear(Matrix& A, vector<U>& b, vector<U>& result);

			U Min(U a, U b);
			U Max(U a, U b);
			U Sign(U a, U b);
			U Sign(U a);
			U Divide(U a, U b);

			ushort_t nDataPoints;
			Data::Data<U> data;

			virtual Function::Function<U>& function() = 0;
			virtual Function::Gradient<U>& gradient() = 0;
			virtual Point::Point<U>& newPoint() = 0;
			virtual Point::Point<U>& deltaPoint() = 0;
			virtual Point::Point<U>& h_sd() = 0;
			virtual Point::Point<U>& h_gn() = 0;
			virtual Point::Point<U>& h_dl() = 0;

			U TINY;
			U zeroLevel;
		};
	}
}
