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
namespace Library
{
	namespace Function
	{
		template <typename T>
		class Function
		{
		public:
			Function(Index N);
			virtual ~Function() { }
			virtual Vector<T>& operator[](Point::Point<T>& other) = 0;
			virtual void updateData() = 0;
			virtual std::string getType() { return "Function.Virtual"; }
		protected:
			std::vector<T> currentData;
			Point::Point<T>* currentPoint;
			Index N;
			bool pointInitialized;
		};

		template <typename T>
		Function<T>::Function(Index N)
		{
			this->N = N;
			for(Index i = 0; i < N; i++)
				currentData.push_back(0.);
			currentPoint = 0;
			pointInitialized = false;
		}
	}
}
