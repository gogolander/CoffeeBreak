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
			virtual Point<T>* clone() = 0;

			virtual string getType() { return string("Point.Virtual"); }

			bool equals(Point<T>& other);
			T& get(ushort_t index) { return point[index]; }
			T& get(string term) { return point[stringToIndex(term)]; }
			string toString();
			ushort_t getDimensions() { return nDimensions; }

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
			if(!other || !this)
				return false;
			if(other.getDimensions() != this->nDimensions)
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
	}
}
