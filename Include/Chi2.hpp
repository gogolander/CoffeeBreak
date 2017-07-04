/*
 * Chi2.hpp
 *
 *  Created on: 20 giu 2017
 *      Author: enzo
 */

#pragma once
#include "Data.hpp"
#include "Matrix.hpp"
#include "Point.hpp"
#include "Maintainance.hpp"
#include <vector>
#include <cstdio>

using namespace std;

namespace Library
{
	namespace Data
	{
		template <typename T>
		class Chi2
		{
		public:
			static T getChi2(Data<T>& data, Point::Point<T>& point, Function::Function<T>& function);
			static vector<T> getGradChi2(Data<T>& data, Point::Point<T>& point,
					Function::Function<T>& function, Function::Gradient<T>& gradient);
			static Matrix<T> getHessianChi2(Data<T>& data, Point::Point<T>& point,
					Function::Function<T>& function, Function::Gradient<T>& gradient);
		};

		template <typename T>
		T Chi2<T>::getChi2(Data<T>& data, Point::Point<T>& point, Function::Function<T>& function)
		{
			if(!point)
				return 1E10;
			T result = 0.;
			for(ushort_t x = 0; x < data.getNDataPoints(); x++)
			{
				T sig2 = 1. / (data.getDeviation()[x] * data.getDeviation()[x]);
				result += sig2 * (data[x] - function[point][x]) *
						(data[x] - function[point][x]);
			}
			return result;
		}

		template <typename T>
		vector<T> Chi2<T>::getGradChi2(Data<T>& data, Point::Point<T>& point,
				Function::Function<T>& function, Function::Gradient<T>& gradient)
		{
			vector<T> gradChi2(point.getDimensions());
			for(ushort_t x = 0; x < data.getNDataPoints(); x++)
			{
				T sig2 = 1. / (data.getDeviation()[x] * data.getDeviation()[x]);
				for(ushort_t i = 0; i < point.getDimensions(); i++)
					gradChi2[i] += (data[x] - function[point][x]) *
						gradient[point][x]->get(i) * sig2;
			}
			return gradChi2;
		}

		template <typename T>
		Matrix<T> Chi2<T>::getHessianChi2(Data<T>& data, Point::Point<T>& point,
				Function::Function<T>& function, Function::Gradient<T>& gradient)
		{
			Matrix<T> hessianChi2(point.getDimensions(), point.getDimensions());
			for(ushort_t x = 0; x < data.getNDataPoints(); x++)
			{
				T sig2 = 1. / (data.getDeviation()[x] * data.getDeviation()[x]);
				for(ushort_t j = 0; j < point.getDimensions(); j++) // Use the simmetry, Luke!
					for(ushort_t i = 0; i <= j; i++)
					{
						T term = sig2 * (gradient[point][x]->get(i) *
								gradient[point][x]->get(j));
						hessianChi2[i][j] += term;
//						printf("[%i,%i]: %G * %G * %G = %G\n", j, i,
//								gradient[point][x]->get(i), gradient[point][x]->get(j), sig2,
//								term);
					}
			}
			for(ushort_t j = 0; j < point.getDimensions(); j++)
				for(ushort_t i = 0; i < j; i++)
					hessianChi2[j][i] = hessianChi2[i][j];
			return hessianChi2;
		}
	}
}
