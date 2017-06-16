/*
 * Function.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Point.hpp"
#include "Maintainance.hpp"
#include <string>

using namespace std;

namespace Library
{
	namespace Function
	{
		template <typename T>
		class Function
		{
		public:
			Function() { currentPoint = 0; nDataPoints = 0; }
			virtual ~Function() { }

			virtual void updateData() = 0;
			virtual string getType() { return "Function.Virtual"; }

			vector<T>& operator[](Point::Point<T>& other);
			T getFlux();
		protected:
			vector<T> currentData;
			Point::Point<T>* currentPoint;
			ushort_t nDataPoints;
		};

		template <typename T>
		vector<T>& Function<T>::operator[](Point::Point<T>& other)
		{
			if(other)
			{
				if(currentPoint != 0 && !currentPoint->equals(other))
				{
					currentPoint = other.clone();
					updateData();
				}
				else if(currentPoint == 0)
				{
					currentPoint = other.clone();
					updateData();
				}
			}
			return currentData;
		}

		template <typename T>
		T Function<T>::getFlux()
		{
			if(currentPoint != 0)
				return 0;
			T result = 0;
			for(ushort_t index = 0; index < nDataPoints; index++)
				result += currentData[index];
			return result;
		}
	}
}
