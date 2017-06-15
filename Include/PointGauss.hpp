/*
 * GaussPoint.hpp
 *  y(x) = a * exp(-1 * (x - c)^2 / (2 * s^2)) + d
 *
 *  FWHM = c * 2 * SQRT(2 * LOG(2))
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */
#pragma once
#include <string>
#include <iostream>

#include "Maintainance.hpp"
#include "Point.hpp"

using namespace std;

namespace Library
{
	namespace Point
	{
		template <typename T>
		class PointGauss : public Point<T>
		{
		public:
			PointGauss()
			{
				this->point = vector<T>(4);
				this->nDimensions = 4;
			}
			PointGauss(PointGauss<T>& other);
			PointGauss(PointGauss<T>* other);

			virtual ~PointGauss() { }
			virtual bool isValid();
			virtual void copyTo(Point<T>* target);
			virtual PointGauss<T> clone();

			virtual PointGauss<T> operator +(const T& b);

			ushort_t stringToIndex(string term);
			string getType() { return "Point.Gauss"; }

			T getFWHM() { return this->get("sigma") / GAUSS_FWHM; }
		};

		template <typename T>
		PointGauss<T>::PointGauss(PointGauss<T>& other) : PointGauss<T>()
		{
			for(ushort_t index = 0; index < this->nDimensions; index++)
				this->get(index) = other[index];
		}

		template <typename T>
		PointGauss<T>::PointGauss(PointGauss<T>* other) : PointGauss<T>()
		{
			for(ushort_t index = 0; index < this->nDimensions; index++)
				this->get(index) = other->get(index);
		}

		template <typename T>
		bool PointGauss<T>::isValid() {	return true; }

		template <typename T>
		PointGauss<T> PointGauss<T>::clone()
		{
			PointGauss<T> result; // I smell trouble
			for(ushort_t index = 0; index < this->nDimensions; index++)
				result[index] = this->get(index);
			return result;
		}

		template <typename T>
		ushort_t PointGauss<T>::stringToIndex(string term)
		{
			transform(term.begin(), term.end(), term.begin(), ::tolower);
			ushort_t value = (term == "amplitude" || term == "amp" || term == "a") * 1 +
					(term == "center" || term == "c") * 2 +
					(term == "s" || term == "sigma" || term == "dispersion" ||
							term == "disp") * 3 +
					(term == "zero" || term == "zerolevel" || term == "zero level" ||
							term == "const") * 4;
			return value - 1;
		}

		template <typename T>
		void PointGauss<T>::copyTo(Point<T>* target)
		{
			if(target == NULL)
				target = new PointGauss<T>();
			for(ushort_t i = 0; i < this->nDimensions; i++)
				target->get(i) = this->point[i];
		}

		template <typename T>
		PointGauss<T> PointGauss<T>::operator +(const T& b)
		{
			PointGauss<T> c = this->clone();
			for(ushort_t index = 0; index < c.getDimensions(); index++)
				c[index] += b;
			return c;
		}
	}
}
