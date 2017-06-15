#pragma once
#include <boost/shared_ptr.hpp>
#include <string>
#include "Point.hpp"
#include "Function.hpp"
#include "Gradient.hpp"

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
			typedef typename boost::shared_ptr<U**> Matrix;

			Fit();
			Fit(Vector<U>& data, int nDataPoints);
			int LevenbergMarquardt(Point::Point<U>& point, U gradientTOL, U xTOL, U chi2TOL, int maxIterations);
			int DogLeg(Point::Point<U>& point, U Delta0, U gradientTOL, U xTOL, U chi2TOL, int maxIterations);
			U Chi2(Point::Point<U>& point); // This really is confusing: simplify it!
			void getData(Vector<U>& result);
			void setData(Vector<U>& data);
			void getModel(Point::Point<U>& point, Vector<U>& result);
			void getResiduals(Point::Point<U>& point, Vector<U>& result);
			std::string PrintPoint(Point::Point<U>& point);
		protected:
			void Chi2Gradient(Point::Point<U>& point, Point::Point<U>& result);
			void Chi2Hessian(Point::Point<U>& point, Matrix& result);
			U Abs(Vector<U>& A, int size);
			void Inverse(Matrix& A, int rank);
			void SolveLinear(Matrix& A, Vector<U>& b, Vector<U>& result);
			U Min(U a, U b);
			U Max(U a, U b);
			U Sign(U a, U b);
			U Sign(U a);
			U Divide(U a, U b);

			int nDataPoints;
			Vector<U> data;
			Vector<U> deviation;
			Vector<U> model;
			Vector<U> residuals;
			Function::Function<U>* function;
			Function::Gradient<U>* gradient;
			U TINY;
			U zeroLevel;
		};
	}
}
