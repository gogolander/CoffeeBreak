/*
 * PointMoffat.hpp
 *	y(x) = a * (1 + b * (x - c)^2) ^ (-beta) + d
 *
 *	FWHM = 2 * sqrt( (2^(1 / beta) - 1) / b )
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include <cmath>
#include "Maintainance.hpp"
#include "Point.hpp"
//#include "PointAGN.hpp"

using namespace std;

namespace Library
{
	namespace Point
	{
		template <typename T>
		class PointMoffat : public Point<T>
		{
		public:
			PointMoffat()
			{
				this->point = vector<T>(5);
				this->nDimensions = 5;
			}
			PointMoffat(PointMoffat<T>& other);
			PointMoffat(PointMoffat<T>* other);

			virtual ~PointMoffat() { }
			virtual bool isValid();
			virtual void copyTo(Point<T>* target);
			virtual PointMoffat<T>* clone() { return new PointMoffat<T>(this); }

			PointMoffat<T> operator +(const T& b);
			PointMoffat<T> operator +(PointMoffat<T>& b);
			PointMoffat<T> operator -(const T& b);
			PointMoffat<T> operator -(PointMoffat<T>& b);
			PointMoffat<T> operator *(const T& b);
			PointMoffat<T> operator *(PointMoffat<T>& b);
			PointMoffat<T> operator /(const T& b);
			PointMoffat<T> operator /(PointMoffat<T>& b);

			ushort_t stringToIndex(string term);
			string getType() { return "Point.Moffat"; }

			T getFWHM() { return 2. * sqrt((pow(2, 1. / this->point["beta"]) - 1) / this->point["b"]); }
		};

		template <typename T>
		PointMoffat<T>::PointMoffat(PointMoffat<T>& other) : PointMoffat<T>()
		{
			for(ushort_t index = 0; index < this->nDimensions; index++)
				this->get(index) = other[index];
		}

		template <typename T>
		bool PointMoffat<T>::isValid()
		{
			for(Index i; i < this->nDimensions; i++)
				if(this->point[i] <= 0)
					return false;
			return true;
		}

		template <typename T>
		ushort_t PointMoffat<T>::stringToIndex(string term)
		{
			transform(term.begin(), term.end(), term.begin(), ::tolower);
			ushort_t value = (term == "amplitude" || term == "amp" || term == "a") * 1 +
					(term == "b") * 2 +
					(term == "center" || term == "c") * 3 +
					(term == "beta") * 4 +
					(term == "zero" || term == "zerolevel" || term == "zero level" ||
							term == "const") * 5;
			return value - 1;
		}

		template <typename T>
		void PointMoffat<T>::copyTo(Point<T>* target)
		{
			if(target == NULL)
				target = new PointMoffat<T>();
			for(ushort_t i = 0; i < this->nDimensions; i++)
				target->get(i) = this->point[i];
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator +(const T& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) + b;
			return c;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator +(PointMoffat<T>& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) + b[i];
			return c;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator -(const T& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) - b;
			return c;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator -(PointMoffat<T>& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) - b[i];
			return c;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator *(const T& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) * b;
			return c;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator *(PointMoffat<T>& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) * b[i];
			return c;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator /(const T& b)
		{
			if(b != 0)
			{
				PointMoffat<T> c;
				for(ushort_t i; i < this->nDimensions; i++)
					c[i] = this->get(i) / b;
				return c;
			}
			return 0;
		}

		template <typename T>
		PointMoffat<T> PointMoffat<T>::operator /(PointMoffat<T>& b)
		{
			PointMoffat<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
			{
				if(b[i] != 0)
					c[i] = this->get(i) / b[i];
				else
					c[i] = -1;
			}
			return c;
		}
	}
}
