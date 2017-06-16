/*
 * PointAGN.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include <string>

#include "Maintainance.hpp"
#include "Point.hpp"

namespace Library
{
	namespace Point
	{
		template <typename T>
		class PointAGN : public Point<T>
		{
		public:
			PointAGN()
			{
				this->point = vector<T>(8);
				this->nDimensions = 8;
			}
			PointAGN(PointAGN<T>& other);
			PointAGN(PointAGN<T>* other);

			virtual ~PointAGN() { }
			virtual bool isValid();
			virtual void copyTo(Point<T>* target);
			virtual PointAGN<T>* clone() { return new PointAGN<T>(this); }

			PointAGN<T> operator +(const T& b);
			PointAGN<T> operator +(PointAGN<T>& b);
			PointAGN<T> operator -(const T& b);
			PointAGN<T> operator -(PointAGN<T>& b);
			PointAGN<T> operator *(const T& b);
			PointAGN<T> operator *(PointAGN<T>& b);
			PointAGN<T> operator /(const T& b);
			PointAGN<T> operator /(PointAGN<T>& b);

			ushort_t stringToIndex(string term);
			string getType() { return "Point.AGN"; }

			T hostFWHM() { return 2. * sqrt((pow(2, 1. / this->point["host beta"]) - 1) / this->point["host b"]); }
			T bhFWHM() { return this->get("bh sigma") / GAUSS_FWHM; }
		};

		template <typename T>
		PointAGN<T>::PointAGN(PointAGN<T>& other) : PointAGN<T>()
		{
			for(ushort_t index = 0; index < this->nDimensions; index++)
				this->get(index) = other[index];
		}

		template <typename T>
		PointAGN<T>::PointAGN(PointAGN<T>* other) : PointAGN<T>()
		{
			for(ushort_t index = 0; index < this->nDimensions; index++)
				this->get(index) = other->get(index);
		}

		template <typename T>
		bool PointAGN<T>::isValid()
		{
			for(ushort_t i; i < this->nDimensions; i++)
				if(this->point[i] <= 0)
					return false;
			return true;
		}

		template <typename T>
		ushort_t PointAGN<T>::stringToIndex(string term)
		{
			transform(term.begin(), term.end(), term.begin(), ::tolower);
			ushort_t value = (term == "bh amplitude" || term == "bh amp" || term == "bh a") * 1 +
					(term == "bh center") * 2 +
					(term == "bh sigma" || term == "bh dispersion") * 3 +
					(term == "host amplitude" || term == "host amp" || term == "host a") * 4 +
					(term == "host b") * 5 +
					(term == "host center" || term == "host c") * 6 +
					(term == "host beta" || term == "host damping") * 7 +
					(term == "zero" || term == "zeroLevel" || term == "zero level" ||
							term == "Zero" || term == "const") * 8;
			return value - 1;
		}

		template <typename T>
		void PointAGN<T>::copyTo(Point<T>* target)
		{
			if(target == NULL)
				target = new PointAGN<T>();
			for(ushort_t i = 0; i < this->nDimensions; i++)
				target->get(i) = this->point[i];
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator +(const T& b)
		{
			PointAGN<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) + b;
			return c;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator +(PointAGN<T>& b)
		{
			PointAGN<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) + b[i];
			return c;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator -(const T& b)
		{
			PointAGN<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) - b;
			return c;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator -(PointAGN<T>& b)
		{
			PointAGN<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) - b[i];
			return c;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator *(const T& b)
		{
			PointAGN<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) * b;
			return c;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator *(PointAGN<T>& b)
		{
			PointAGN<T> c;
			for(ushort_t i; i < this->nDimensions; i++)
				c[i] = this->get(i) * b[i];
			return c;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator /(const T& b)
		{
			if(b != 0)
			{
				PointAGN<T> c;
				for(ushort_t i; i < this->nDimensions; i++)
					c[i] = this->get(i) / b;
				return c;
			}
			return 0;
		}

		template <typename T>
		PointAGN<T> PointAGN<T>::operator /(PointAGN<T>& b)
		{
			PointAGN<T> c;
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
