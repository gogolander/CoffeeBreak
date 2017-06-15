/*
 * Point.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include <string>
#include <vector>
#include <iostream>

#include "Maintainance.hpp"

using namespace std;

namespace Library
{
	namespace Point
	{
		template <typename T>
		class Point
		{
		public:
			virtual ~Point() { }
			virtual bool isValid() = 0;
			virtual ushort_t stringToIndex(string term) = 0;
			virtual void copyTo(Point<T>* target) = 0;
			virtual string getType() { return "Point.Virtual"; }
			virtual Point<T> clone() = 0;

			bool equals(Point<T>& other);
			T& get(ushort_t index) { return point[index]; }
			T& get(string term) { return point[stringToIndex(term)]; }
			string toString();
			ushort_t getDimensions() { return nDimensions; }

			Point<T> operator +(const Point<T>& b);
			virtual Point<T> operator +(const T& b) = 0;

			bool operator ==(const Point<T>& other) { return this->equals(other); }
			bool operator !=(const Point<T>& other) { return !this->equals(other); }

			T& operator [](ushort_t param) { return get(param); }
			T& operator [](const char param[]) { return get(string(param)); }
			T& operator [](string term) { return get(term); }

			operator string()
			{
				string result = "";
				for(ushort_t param = 0; param < nDimensions; param++)
					result += to_string(point[param]) + " ";
				return result;
			}

			operator bool()
			{
				if(this == 0)
					return false;
				return isValid();
			}
		protected:
			vector<T> point;
			ushort_t nDimensions;
		};

		template <typename T>
		bool Point<T>::equals(Point<T>& other)
		{
			if(!other)
				return false;
			if(!this)
				return false;
			if(other.getType() != this->getType())
				return false;
			for(ushort_t param = 0; param < this->nDimensions; param++)
				if(other[param] != this->point[param])
					return false;
			return true;
		}

		template <typename T>
		string Point<T>::toString()
		{
			string result = "";
			for(ushort_t index = 0; index < nDimensions; index++)
				result += std::to_string(point[index]) + " ";
			return result;
		}

		template <typename T>
		Point<T> Point<T>::operator +(const Point<T>& b)
		{
			Point<T> c = b;
			for(ushort_t index; index < c.getDimensions(); index++)
				b[index] += get(index);
			return c;
		}
	}
}
