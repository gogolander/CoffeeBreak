#include "Point.hpp"
#include "Function.hpp"
#include "Gradient.hpp"
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include "Fit.h"
#include "Maintainance.hpp"
#include "Data.hpp"
#include "Chi2.hpp"
#include "Matrix.hpp"

using namespace std;

namespace Library
{
	namespace Fit
	{
		/*
		 * Costruttori
		 */
		template <typename U>
		Fit<U>::Fit()
		{
			TINY = 1E-10;
			nDataPoints = 0;
		}

		/*
		 * Getters and Setters
		 */
		template <typename U>
		Data::Data<U>& Fit<U>::getData() { return this->data; }

		template <typename U>
		void Fit<U>::setData(vector<U>& data)
		{
			this->data = Data::Data<U>(data, data.size());
			this->data.initData();
		}

		template <typename U>
		vector<U>& Fit<U>::getModel(Point::Point<U>& point)	{ return this->function[point]; }

		template <typename U>
		vector<U>& Fit<U>::getResiduals(Point::Point<U>& point)
		{
			vector<U> result;
			for(ushort_t i = 0; i < nDataPoints; i++)
				result[i] = this->data[i] - this->function[point][i];
			return result;
		}
		/*
		 * General Purpose
		 */
		template <typename U>
		U Fit<U>::Min(U a, U b)
		{
			if(a > b) return b;
			return a;
		}

		template <typename U>
		U Fit<U>::Max(U a, U b)
		{
			if(a > b) return a;
			return b;
		}

		template <typename U>
		U Fit<U>::Sign(U a, U b) { return Sign(a) * fabs(b); }

		template <typename U>
		U Fit<U>::Sign(U a)
		{
			if(a < 0) return -1;
			return 1;
		}

		template <typename U>
		U Fit<U>::Divide(U a, U b)
		{
			if(a == 0 && b != 0) return 0;
			if(a != 0 && b == 0) return Sign(a) * 1E6;
			if(a == 0 && b == 0) return 0;
			if(b == 1) return a;
			if(b == -1) return -1 * a;
			return a / b;
		}

		template <typename U>
		U Fit<U>::abs(vector<U>& vec)
		{
			U result = 0;
			for(ushort_t i = 0; i < vec.size(); i++)
				result += vec[i] * vec[i];
			return result;
		}

//		template <typename U>
//		void Fit<U>::SolveLinear(Matrix& A, vector<U>& b, vector<U>& result)
//		{
//			ushort_t nParameters = sizeof(b) / sizeof(U);
//			int* index = new int[nParameters];
//			Matrix B(boost::extents[nParameters][nParameters + 1]);
//			typedef typename Matrix::index dimensions;
//
//			// *** Inizializzo la matrice dei fattori *** //
//			for(dimensions i = 0; i < nParameters; i++)
//			{
//				index[i] = i;
//				for(dimensions j = 0; j < nParameters; j++)
//					B[i][j] = A[i][j];
//				B[i][nParameters] = b[i];
//			}
//			for(dimensions i = 0; i < nParameters; i++)
//			{
//				// *** Remove zeros from the principal diagonal *** //
//				if(B[i][i] == 0)
//				{
//					// *** Search a row without zero on the diagonal *** //
//					for(dimensions j = 0; j < nParameters; j++)
//					{
//						if(B[i][j] != 0 )
//						{
//							// *** Swap the i-th row with the j-th row *** //
//							for(dimensions k = 0; k <= nParameters; k++)
//							{
//								U swap = B[k][i];
//								B[k][i] = B[k][j];
//								B[k][j] = swap;
//								// *** remember to swap the map to retrieve solutions too! *** //
//								swap = index[i];
//								index[i] = index[j];
//								index[j] = (int)swap;
//							}
//							// *** Goal accomplished *** //
//							break;
//						}
//					}
//				}
//				// *** Pivoting routine *** //
//				if(B[i][i] != 1)
//					for(dimensions k = 0; k <= nParameters; k++)
//						if(B[i][k] != 0) // Skip the zeros
//							B[i][k] = Divide(B[i][k], B[i][i]);
//				// *** Eliminate elements below the pivots *** //
//				U factor = B[i][i + 1];
//				for(dimensions j = i + 1; j < nParameters; j++)
//				{
//					U factor = B[i][j];
//					for(dimensions k = 0; k <= nParameters; k++)
//						B[j][k] -= factor * B[i][k];
//				}
//			}
//			// *** Eliminate elements above the pivot *** //
//			for(dimensions j = nParameters - 1; j > 0; j--)
//			{
//				for(dimensions i = j - 1; i >= 0; i--)
//				{
//					U factor = B[i][j];
//					for(dimensions k = j; k <= nParameters; k++)
//						if(B[j][k] != 0)
//							B[i][k] -= B[j][k] * factor;
//				}
//			}
//			// *** Copy the solution *** //
//			for(dimensions i = 0; i < nParameters; i++)
//				result[i] = B[index[i]][nParameters];
//		}

		/*
		 * Fit to a model
		 */
		/// <summary><para>Fit using the Levenberg-Marquardt algorithm.</para>
		///			<para>The data must were already given to the constructor.</para></summary>
		/// <param name="Point">Initial guess. Contains the final result of the fit.</param>
		/// <param name="gradientTOL">Threshold below which consider the gradient negligible. Default: 1E-8</param>
		/// <param name="xTOL">Threshold below which consider the point shift negligeble and the point stable. Default: 1E-8</param>
		/// <param name="chiTOL">Threshold below which consider the fit acceptable. Default: 1E-4</param>
		/// <param name="maxIterations">Max iterations allowed. Default: 1000</param>
		/// <returns>-1: fail; 0: chi2 negligeble; 1: gradient negligeble; 2: point shift negligeble</returns>
		template <typename U>
		int Fit<U>::LevenbergMarquardt(Point::Point<U>& point, U gradientTOL, U xTOL, U chi2TOL, int maxIterations)
		{
//			for(ushort_t index = 0; index < nDataPoints; index++)
//				std::cout << this->data[index] << " " << this->deviation[index] << std::endl;
			U l = 0;
			U oldChi2 = 0., newChi2 = 0.;
			ushort_t range = point.getDimensions();
			newPoint() = point;
			oldChi2 = Data::Chi2<U>::getChi2(data, point, function());
			vector<U> gradChi2(Data::Chi2<U>::getGradChi2(data, point, function(), gradient())); // Se ci metto il gradiente della funzione non fa altro che trovare il vertice
			Data::Matrix<U> hessianChi2;
			hessianChi2 = Data::Chi2<U>::getHessianChi2(data, point, function(), gradient());
			l = 0.001;
			for(int iter = 0; iter < maxIterations; iter++)
			{
				if(abs(gradChi2) < gradientTOL)
				{
					point = newPoint();
					cout << "Gradient small enough.\nBest fit found in " << iter << " iterations" << endl;
					return 1;
				}

				// *** Determine new LM step *** //
				for(ushort_t i = 0; i < point.getDimensions(); i++)
					hessianChi2[i][i] *= 1.0 + l;

				deltaPoint().fromVector(hessianChi2.solveLinear(gradChi2));
				for(ushort_t j = 0; j < point.getDimensions(); j++)
					newPoint()[j] = point[j] + deltaPoint()[j];
				if(deltaPoint().abs() < xTOL)
				{
					cout << "Shift small enough.\nBest fit found in " << iter << " iterations" << endl;
					point = newPoint();
					return 2;
				}
				newChi2 = Data::Chi2<U>::getChi2(data, newPoint(), function());
//				cout << "Point: " << (string)point << endl;
//				cout << "New Point: " << (string)newPoint() << endl;
//				cout << "Delta Point: " << (string)deltaPoint() << endl;
//				printf("lambda: %G\toldChi2: %G\tnewChi2: %G\n", l, oldChi2, newChi2);
				if(oldChi2 > newChi2)
				{
					point = newPoint();
					if(fabs(oldChi2 - newChi2) < chi2TOL)
					{
						cout << "Chi2 variation small enough.\nBest fit found in " << iter << " iterations." << endl;
						return 0;
					}
					oldChi2 = newChi2;
					gradChi2 = Data::Chi2<U>::getGradChi2(data, point, function(), gradient());
					hessianChi2 = Data::Chi2<U>::getHessianChi2(data, point, function(), gradient());
					l *= 0.1;
				}
				else
				{
					if(2 * fabs(l) > numeric_limits<U>::max())
					{
						cout << "Overflow on mu in the next iteration: " << l << " (l: " << numeric_limits<double>::max() << ")\n";
						break;
					}
					l *= 5; // perturba l ma resta in zona, non andare troppo verso l'ignoto
				}
			}
			cout << "Fit not found" << endl;
			return -1;
		}

		template <typename U>
		int Fit<U>::DogLeg(Point::Point<U>& point, U Delta0, U gradientTOL, U xTOL, U chi2TOL, int maxIterations)
		{
			U alpha = 0, Delta = Delta0, beta = 0, rho, newChi2, oldChi2;
			oldChi2 = Data::Chi2<U>::getChi2(data, point, function());
			vector<U> gradChi2;
			Data::Matrix<U> hessianChi2;
			for(int iter = 0; iter < maxIterations; iter++)
			{
				gradChi2(Data::Chi2<U>::getGradChi2(data, point, function(), gradient()));
				hessianChi2 = Data::Chi2<U>::getHessianChi2(data, point, function(), gradient());
				hessianChi2.inverse();
				for(ushort_t i = 0; i < point.getDimensions(); i++)
				{
					h_gn()[i] = 0;
					for(ushort_t j = 0; j < point.getDimensions(); j++)
						h_gn()[i] -= hessianChi2[j][i] * gradChi2()[j];
					h_sd()[i] = -1 * gradChi2()[i];
				}
				alpha = 0;
				for(ushort_t i = 0; i < point.getDimensions(); i++)
					for(int x = 0; x < nDataPoints; x++)
						alpha += gradient()[point][x][i] * gradient()[point][x][i] *
							gradChi2()[i] * gradChi2()[i];

				alpha = gradChi2().abs() / alpha;
				if(h_gn().abs() <= Delta)
				{
					for(int i = 0; i < point.getDimensions(); i++)
						h_dl()[i] = h_gn()[i];
					rho = oldChi2;
				}
				else if(alpha * h_sd().abs() >= Delta)
				{
					U abs = h_sd().abs();
					for(int i = 0; i < point.getDimensions(); i++)
						h_dl()[i] = Delta * h_sd()[i] / abs;
					rho = Delta * (2 * alpha * gradChi2().abs() - Delta) / (2 * alpha);
				}
				else
				{
					U c = 0, ba = 0, a = 0;
					for(int i = 0; i < point.getDimensions(); i++)
					{
						a += (alpha * alpha * h_sd()[i] * h_sd()[i]);
						ba += (h_gn()[i] - alpha * h_sd()[i]) * (h_gn()[i] - alpha * h_sd()[i]);
						c += alpha * h_gn()[i] * (h_sd()[i] - alpha * h_gn()[i]);
					}
					if(c <= 0)
						beta = (-1 * c + sqrt(c * c + ba * (Delta * Delta - a))) / ba;
					else
						beta = (Delta * Delta - a) / (c + sqrt(c * c + ba *(Delta * Delta - a)));
					for(int i = 0; i < point.getDimensions(); i++)
						h_dl()[i] = alpha * h_sd()[i] + beta * (h_gn()[i] - alpha * h_sd()[i]);
					U tmp = gradChi2().abs();
					rho = 0.5 * alpha * (1 - beta) * (1 - beta) * tmp * tmp +
							beta * (2 - beta) * oldChi2;
				}
				if(h_dl().abs() < xTOL * (point.abs() + xTOL))
				{
					cout << "Best fit found in " << iter << " iterations." << endl;
					return 1;
				}
				for(int i = 0; i < point.getDimensions(); i++)
					newPoint()[i] = point[i] + h_dl()[i];
				newChi2 = Data::Chi2<U>::getChi2(data, newPoint(), function());
				rho = (oldChi2 - newChi2) / rho;
				if(rho > 0 && newChi2 < oldChi2)
				{
					for(int i = 0; i < point.getDimensions(); i++)
						point[i] = newPoint()[i];
					gradChi2 = Data::Chi2<U>::getGradChi2(data, point, function(), gradient());
					hessianChi2 = Data::Chi2<U>::getHessianChi2(data, point, function(), gradient());
					oldChi2 = newChi2;
					if(oldChi2 < chi2TOL)
					{
						cout << "Best fit found in " << iter << " iterations." << endl;
						return 0;
					}
					if(gradChi2().abs() < gradientTOL)
					{
						cout << "Best fit found in " << iter << " iterations." << endl;
						return 2;
					}
				}
				if(rho > 0.75)
				{
					double tmp = 3 * h_dl().abs();
					if(Delta < tmp)
						Delta = tmp;
				}
				else if(rho < 0.25)
				{
					Delta *= 0.5;
					if(Delta < xTOL * (point.abs() + xTOL))
					{
						cout << "Best fit found in " << iter << " iterations." << endl;
						return 2;
					}
				}
			}
			cout << "Fit not found" << endl;
			return -1;
		}
	}
}
