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
namespace Library
{
	namespace Function
	{
		template <typename T>
		class Gradient
		{
		public:
			Gradient(Index N);
			virtual ~Gradient() { }
			virtual Vector<Point::Point<T>*>& operator[](Point::Point<T>& other) = 0;
			virtual void updateData() = 0;
			virtual std::string getType() { return "Gradient.Virtual"; }
		protected:
			std::vector<Point::Point<T>*> currentData;
			Point::Point<T>* currentPoint;
			Index N;
			bool pointInitialized;
		};

		template <typename T>
		Gradient<T>::Gradient(Index N)
		{
			this->N = N;
			for(Index i = 0; i < N; i++)
				currentData.push_back(0.);
			currentPoint = 0;
			pointInitialized = false;
		}
	}
}
