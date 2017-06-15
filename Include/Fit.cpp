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
#include <boost/multi_array.hpp>
#include "Fit.h"
#include "Maintainance.hpp"

using namespace std;

namespace Library
{
	namespace Fit
	{
		/*
		 * Costruttori
		 */
		template <typename U>
		Fit<U>::Fit() { TINY = 1E-10; nDataPoints = 0; }

		template <typename U>
		Fit<U>::Fit(Vector<U>& data, int nDataPoints)
		{
			this->TINY = 1E-10;
			this->nDataPoints = nDataPoints;
			this->data(new U[nDataPoints]);
			this->deviation(new U[nDataPoints]);
			this->model(new U[nDataPoints]);
			this->residuals(new U[nDataPoints]);
			U max = data[0], y = 0, y2 = 0, points = 0, gain = 0.68;
			int semi_width = 12;
			U center;
			for(int i = 0; i < nDataPoints; i++)
			{
				if(data[i] > max)
				{
					max = data[i];
					center = i;
				}
			}
			for(int i = 0; i < nDataPoints; i++)
			{
				if(i < center - semi_width || i > center + semi_width)
				{
					y += data[i];
					y2 += data[i] * data[i];
					points++;
				}
			}
			U var = (y2 / points) - y * y / (points * points);
			for(int i = 0; i < nDataPoints; i++)
			{
				if(i > center - semi_width && i < center + semi_width)
				{
					this->data[i] = data[i];
					this->deviation[i] = sqrt(fabs(data[i]) / gain + var);
				}
				else
				{
					this->data[i] = sqrt(var);
					this->deviation[i] = sqrt(var);
				}
			}
			this->zeroLevel = sqrt(var);
		}

		/*
		 * Getters and Setters
		 */
		template <typename U>
		void Fit<U>::getData(Vector<U>& result)
		{
			for(Index i = 0; i < nDataPoints; i++)
				result[i] = this->data[i];
		}
		template <typename U>
		void Fit<U>::setData(Vector<U>& data)
		{
			for(Index i = 0; i < nDataPoints; i++)
				this->data[i] = data[i];
		}
		template <typename U>
		void Fit<U>::getModel(Point::Point<U>& point, Vector<U>& result)
		{
			for(Index i = 0; i < nDataPoints; i++)
				result[i] = this->function[point][i];
		}
		template <typename U>
		void Fit<U>::getResiduals(Point::Point<U>& point, Vector<U>& result)
		{
			for(Index i = 0; i < nDataPoints; i++)
				result[i] = this->data[i] - this->function[point][i];
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
		U Fit<U>::Sign(U a, U b)
		{
			return Sign(a) * fabs(b);
		}
		template <typename U>
		U Fit<U>::Sign(U a)
		{
			if(a < 0)
				return -1;
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
		string Fit<U>::PrintPoint(Point::Point<U>& Point)
		{
			std::ostringstream strs;
			strs.precision(5);
			for(int i = 0; i < Point.getDimensions(); i++)
				strs << Point[i] << " ";
			return strs.str();
		}
		/*
		 * Basic matrix algebra
		 */
		template <typename U>
		U Fit<U>::Abs(Vector<U>& A, int size)
		{
			U result = 0.;
			for(int i = 0; i < size; i++)
				result += A[i] * A[i];
			return sqrt(result);
		}
		template <typename U>
		void Fit<U>::Inverse(Matrix& A, int rank)
		{
			// *** Inizializzazione della matrice B = [A|I] *** //
			Matrix B(boost::extents[rank][2 * rank]);
			typedef typename Matrix::index index;

			for(index i = 0; i < rank; i++)
			{
				for(index j = 0; j < rank; j++)					// *** Popolo la metà sinistra di B con la matrice A *** //
					B[i][j] = A[i][j];
				for(index j = rank; j < 2 * rank; j++)
					B[i][j] = (int)(j - rank == i);				// *** Popolo la metà destra di B con la matrice identita' *** //
			}
			// *** Azzero tutti gli elementi fuori diagonale *** //
			for(index k = 0; k < rank; k++)
			{
				if(B[k][k] == 0)								// *** L'elemento sulla diagonale non può mai essere zero! *** //
				{
					for(index i = k + 1; i < rank; i++)
						if(B[i][k] != 0)						// *** Scambia la i-esima riga con la k-esima *** //
						{
							for(index j = 0; j < 2 * rank; j++)
							{
								U swap = B[i][j];
								B[i][j] = B[k][j];
								B[k][j] = swap;
							}
							break;								// *** Fatto. Passa al calcolo dell'inversa *** //
						}
				}
				if(B[k][k] != 1)								// *** Divido in modo da avere l'elemento sulla diagonale pari a 1*** //
					for(index j = 2 * rank - 1; j >= k; j--)
						B[k][j] = Divide(B[k][j], B[k][k]);		// *** Divisione stabilizzata *** //
				for(index i = 0; i < rank; i++)					// *** Sottraggo la k-esima riga alla i-esima *** //
					if(k != i)
						for(index j = 2 * rank - 1; j >= k; j--)	// *** Restringo agli elementi che stanno sulla diagonale *** //
							B[i][j] -= B[i][k] * B[k][j];		// *** Quelli a sx di essa sono già nulli a questo punto *** //
			}
			// *** Copio l'inversa in A *** //
			for(index i = 0; i < rank; i++)
				for(index j = rank; j < 2 * rank; j++)
					A[i][j - rank] = B[i][j];
		}
		template <typename U>
		void Fit<U>::SolveLinear(Matrix& A, Vector<U>& b, Vector<U>& result)
		{
			ushort nParameters = sizeof(b) / sizeof(U);
			int* index = new int[nParameters];
			Matrix B(boost::extents[nParameters][nParameters + 1]);
			typedef typename Matrix::index dimensions;

			// *** Inizializzo la matrice dei fattori *** //
			for(dimensions i = 0; i < nParameters; i++)
			{
				index[i] = i;
				for(dimensions j = 0; j < nParameters; j++)
					B[i][j] = A[i][j];
				B[i][nParameters] = b[i];
			}
			for(dimensions i = 0; i < nParameters; i++)
			{
				// *** Remove zeros from the principal diagonal *** //
				if(B[i][i] == 0)
				{
					// *** Search a row without zero on the diagonal *** //
					for(dimensions j = 0; j < nParameters; j++)
					{
						if(B[i][j] != 0 )
						{
							// *** Swap the i-th row with the j-th row *** //
							for(dimensions k = 0; k <= nParameters; k++)
							{
								U swap = B[k][i];
								B[k][i] = B[k][j];
								B[k][j] = swap;
								// *** remember to swap the map to retrieve solutions too! *** //
								swap = index[i];
								index[i] = index[j];
								index[j] = (int)swap;
							}
							// *** Goal accomplished *** //
							break;
						}
					}
				}
				// *** Pivoting routine *** //
				if(B[i][i] != 1)
					for(dimensions k = 0; k <= nParameters; k++)
						if(B[i][k] != 0) // Skip the zeros
							B[i][k] = Divide(B[i][k], B[i][i]);
				// *** Eliminate elements below the pivots *** //
				U factor = B[i][i + 1];
				for(dimensions j = i + 1; j < nParameters; j++)
				{
					U factor = B[i][j];
					for(dimensions k = 0; k <= nParameters; k++)
						B[j][k] -= factor * B[i][k];
				}
			}
			// *** Eliminate elements above the pivot *** //
			for(dimensions j = nParameters - 1; j > 0; j--)
			{
				for(dimensions i = j - 1; i >= 0; i--)
				{
					U factor = B[i][j];
					for(dimensions k = j; k <= nParameters; k++)
						if(B[j][k] != 0)
							B[i][k] -= B[j][k] * factor;
				}
			}
			// *** Copy the solution *** //
			for(dimensions i = 0; i < nParameters; i++)
				result[i] = B[index[i]][nParameters];
		}
		/*
		 * Chi2
		 */
		template <typename U>
		U Fit<U>::Chi2(Point::Point<U>& point)
		{
			if(!point.isValid())
				return 1E10;
			U result = 0.;
			for(Index x = 0; x < this->nDataPoints; x++)
			{
				U sig2 = 1. / (this->deviation[x] * this->deviation[x]);
				result += sig2 * (this->data[x] - this->function[point][x]) * (this->data[x] - this->function[point][x]);
			}
			return result;
		}
		template <typename U>
		void Fit<U>::Chi2Gradient(Point::Point<U>& point, Point::Point<U>& result)
		{
			for(Index i = 0; i < result.getDimensions(); i++)
				result[i] = 0.;
			for(Index x = 0; x < nDataPoints; x++)
			{
				U sig2 = 1. / (deviation[x] * deviation[x]);
				for(Index i = 0; i < result.getDimensions(); i++)
					result[i] += (data[x] - this->function[point][x]) * this->gradient[x][i] * sig2;
			}
		}
		template <typename U>
		void Fit<U>::Chi2Hessian(Point::Point<U>& point, Matrix& result)
		{
			for(Index j = 0; j < point.getDimensions(); j++) // Use the simmetry, Luke!
				for(Index i = 0; i <= j; i++)
					result[i][j] = 0.;
			for(Index x = 0; x < nDataPoints; x++)
			{
				U sig2 = 1. / (deviation[x] * deviation[x]);
				for(Index j = 0; j < point.getDimensions(); j++) // Use the simmetry, Luke!
					for(Index i = 0; i <= j; i++)
						result[i][j] += sig2 * (this->gradient[x][i] * this->gradient[x][j]);
			}
			for(Index j = 0; j < point.getDimensions(); j++)
				for(Index i = 0; i < j; i++)
					result[j][i] = result[i][j];
		}
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
			U l = 0;
			U oldChi2 = 0., newChi2 = 0.;
			Point::Point<U> gradChi2(point.getDimensions());
			Matrix hessianChi2(boost::extents[point.getDimensions()][point.getDimensions()]);
			typedef typename Matrix::index index;

			boost::shared_ptr<Point::Point<U> > deltaPoint;
			boost::shared_ptr<Point::Point<U> > newPoint;
			point.newInstance(deltaPoint);
			point.newInstance(newPoint);

			for(index i = 0; i < point.getDimensions(); i++)
			{
				deltaPoint[i] = 0.;
				newPoint[i] = point[i];
			}
			oldChi2 = Chi2(point);
			Chi2Gradient(point, gradChi2); // Se ci metto il gradiente della funzione non fa altro che trovare il vertice
			Chi2Hessian(point, hessianChi2);
			l = 0.001;
			for(int iter = 0; iter < maxIterations; iter++)
			{
				if(Abs(gradChi2, point.getDimensions()) < gradientTOL)
				{
					for(int i = 0; i < point.getDimensions(); i++)
						point[i] = newPoint[i];
					cout << "Gradient small enough.\nBest fit found in " << iter << " iterations" << endl;
					return 1;
				}
				for(index i = 0; i < point.getDimensions(); i++)
					hessianChi2[i][i] *= 1.0 + l;
				Inverse(hessianChi2, point.getDimensions());
				for(index i = 0; i < point.getDimensions(); i++)
				{
					deltaPoint[i] = 0;
					for(index j = 0; j < point.getDimensions(); j++)
						deltaPoint[i] += hessianChi2[i][j] * gradChi2[j];
					for(index j = 0; j < point.getDimensions(); j++)
						newPoint[j] = point[j] + deltaPoint[j];
				}
				if(Abs(deltaPoint, point.getDimensions()) < xTOL)
				{
					cout << "Shift small enough.\nBest fit found in " << iter << " iterations" << endl;
					for(index i = 0; i < point.getDimensions(); i++)
						point[i] = newPoint[i];
					return 2;
				}
				newChi2 = Chi2(newPoint);
				if(oldChi2 > newChi2)
				{
					for(index i = 0; i < point.getDimensions(); i++)
						point[i] = newPoint[i];
					if(fabs(oldChi2 - newChi2) < chi2TOL)
					{
						cout << "Chi2 variation small enough.\nBest fit found in " << iter << " iterations." << endl;
						return 0;
					}
					oldChi2 = newChi2;
					Chi2Gradient(point, gradChi2);
					Chi2Hessian(point, hessianChi2);
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
			Matrix hessianChi2(boost::extents[point.getDimensions()][point.getDimensions()]);
			typedef typename Matrix::index index;
			Point::Point<U> gradChi2(point.getDimensions());
			Point::Point<U> h_sd(point.getDimensions()); // Gauss-Newton step
			Point::Point<U> h_gn(point.getDimensions()); // Steepest Descent step
			Point::Point<U> h_dl(point.getDimensions()); // Dog Leg step, aka thed actual step to use
			boost::shared_ptr<Point::Point<U> > newPoint = point.getInstance();
			U alpha = 0, Delta = Delta0, beta = 0, rho, newChi2, oldChi2;
			oldChi2 = this->Chi2(point);
			for(int iter = 0; iter < maxIterations; iter++)
			{
				this->Chi2Gradient(point, gradChi2);
				this->Inverse(hessianChi2, point.getDimensions());
				for(index i = 0; i < point.getDimensions(); i++)
				{
					h_gn[i] = 0;
					for(index j = 0; j < point.getDimensions(); j++)
						h_gn[i] -= hessianChi2[j][i] * gradChi2[j];
					h_sd[i] = -1 * gradChi2[i];
				}
				alpha = 0;
				for(index i = 0; i < point.getDimensions(); i++)
					for(int x = 0; x < nDataPoints; x++)
						alpha += gradient[point][x][i] * gradient[point][x][i] * gradChi2[i] * gradChi2[i];

				alpha = this->Abs(gradChi2, point.getDimensions()) / alpha;
				if(this->Abs(h_gn, point.getDimensions()) <= Delta)
				{
					for(int i = 0; i < point.getDimensions(); i++)
						h_dl[i] = h_gn[i];
					rho = oldChi2;
				}
				else if(alpha * this->Abs(h_sd, point.getDimensions()) >= Delta)
				{
					U abs = this->Abs(h_sd, point.getDimensions());
					for(int i = 0; i < point.getDimensions(); i++)
						h_dl[i] = Delta * h_sd[i] / abs;
					rho = Delta * (2 * alpha * this->Abs(gradChi2, point.getDimensions()) - Delta) / (2 * alpha);
				}
				else
				{
					U c = 0, ba = 0, a = 0;
					for(int i = 0; i < point.getDimensions(); i++)
					{
						a += (alpha * alpha * h_sd[i] * h_sd[i]);
						ba += (h_gn[i] - alpha * h_sd[i]) * (h_gn[i] - alpha * h_sd[i]);
						c += alpha * h_gn[i] * (h_sd[i] - alpha * h_gn[i]);
					}
					if(c <= 0)
						beta = (-1 * c + sqrt(c * c + ba * (Delta * Delta - a))) / ba;
					else
						beta = (Delta * Delta - a) / (c + sqrt(c * c + ba *(Delta * Delta - a)));
					for(int i = 0; i < point.getDimensions(); i++)
						h_dl[i] = alpha * h_sd[i] + beta * (h_gn[i] - alpha * h_sd[i]);
					U tmp = this->Abs(gradChi2, point.getDimensions());
					rho = 0.5 * alpha * (1 - beta) * (1 - beta) * tmp * tmp + beta * (2 - beta) * oldChi2;
				}
				if(this->Abs(h_dl, point.getDimensions()) < xTOL *
						(this->Abs(point, point.getDimensions()) + xTOL))
				{
					cout << "Best fit found in " << iter << " iterations." << endl;
					return 1;
				}
				for(int i = 0; i < point.getDimensions(); i++)
					newPoint[i] = point[i] + h_dl[i];
				newChi2 = this->Chi2(newPoint);
				rho = (oldChi2 - newChi2) / rho;
				if(rho > 0 && newChi2 < oldChi2)
				{
					for(int i = 0; i < point.getDimensions(); i++)
						point[i] = newPoint[i];
					this->Chi2Gradient(point, gradChi2);
					this->Chi2Hessian(point, hessianChi2);
					oldChi2 = newChi2;
					if(oldChi2 < chi2TOL)
					{
						cout << "Best fit found in " << iter << " iterations." << endl;
						return 0;
					}
					if(this->Abs(gradChi2, point.getDimensions()) < gradientTOL)
					{
						cout << "Best fit found in " << iter << " iterations." << endl;
						return 2;
					}
				}
				if(rho > 0.75)
				{
					double tmp = 3 * this->Abs(h_dl, point.getDimensions());
					if(Delta < tmp)
						Delta = tmp;
				}
				else if(rho < 0.25)
				{
					Delta *= 0.5;
					if(Delta < xTOL * (this->Abs(point, point.getDimensions()) + xTOL))
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
