/*
 * Gradient.hpp
 *
 *  Created on: 25 mag 2017
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
		class Gradient
		{
		public:
			Gradient();

			virtual ~Gradient() { }
			virtual void updateData() = 0;
			virtual string getType() { return string("Gradient.Virtual"); }

			vector<Point::Point<T>* >& get(Point::Point<T>& other);
			vector<Point::Point<T>* >& operator[](Point::Point<T>& other) { return get(other); }
		protected:
			vector<Point::Point<T>* > currentGradient;
			Point::Point<T>* currentPoint;
			ushort_t nDataPoints;
		};

		template <typename T>
		Gradient<T>::Gradient()
		{
			nDataPoints = 0;
			currentPoint = 0;
		}

		template <typename T>
		vector<Point::Point<T>* >& Gradient<T>::get(Point::Point<T>& other)
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
			return currentGradient;
		}
	}
}
