#pragma once
#include <complex>
#include <boost/shared_ptr.hpp>
#include "Maintainance.hpp"

using namespace std;
namespace Library
{
	namespace FourierAnalysis
	{
		template <typename T>
		class Fourier
		{
		public:
			Fourier(Vector<T>& data, ushort N);
			void getData(Vector<T>& result);
			void getFftData(Vector<T>& realPart, Vector<T>& imagPart = NULL);
			void getFftData(Vector<complex<T> >& result);
			void getIfftData(Vector<T>& realPart, Vector<T>& imagPart = NULL);
			void getIfftData(Vector<complex<T> >& result);
			ushort log2(ushort N);
			ushort check(ushort n) { return (n > 0 && (n & (n - 1)) == 0); }
			ushort reverse(ushort N, ushort n);
			void ConvolveWith(Vector<T>& filter, ushort N,
					Vector<T>& realResult, Vector<T>& imaginaryResult);
			void ConvolveWith(Fourier<T>& filter, ushort N, Vector<T>& realResult, Vector<T>& imaginaryResult);

			static void FFT(Vector<T>& realPart, Vector<T>& imagPart, ushort N,
					T normalizationFactor);
			static void FFT(Vector<complex<T> >& data, ushort N, T normalizationFactor);
			static void IFFT(Vector<T>& realPart, Vector<T>& imagPart, ushort N);
			static void IFFT(Vector<complex<T> >& data, ushort N);
			static void Convolve(Vector<T>& data, Vector<T>& filter, ushort N,
					Vector<T>& realResult, Vector<T>& imaginaryResult);
		private:
			Vector<T> data;
			Vector<complex<T> > fftData;
			Vector<complex<T> > ifftData; // dovrebbero essere reali ma non si sa mai
			ushort N;
		};

		template <typename T>
		Fourier<T>::Fourier(Vector<T>& data, ushort N)
		{
			this->N = N;
			fftData(N);
			ifftData(N);
			for(Index i = 0; i < N; i++)
				this->fftData[i] = this->data[i] = data[i];
			Fourier<T>::FFT(fftData, N, 1. / N);
			for(Index i = 0; i < N; i++)
				this->ifftData[i] = this->fftData[i];
			Fourier<T>::IFFT(ifftData, N);
		}

		// funzione per calcolare il logaritmo in base 2 di un intero
		template <typename T>
		ushort Fourier<T>::log2(ushort N)
		{
			int k = N, i = 0;
			while(k)
			{
				k >>= 1;
				i++;
			}
			return i - 1;
		}

		// calcola il reverse number di ogni intero n rispetto al numero massimo N
		template <typename T>
		ushort Fourier<T>::reverse(ushort N, ushort n)
		{
			int p = 0;
			int log2N = log2(N);
			for(Index j = 1; j <= log2N; j++)
				if(n & (1 << (log2N - j)))
					p |= 1 << (j - 1);
			return p;
		}

		template <typename T>
		void Fourier<T>::ConvolveWith(Vector<T>& filterData, ushort N,
				Vector<T>& realResult, Vector<T>& imaginaryResult)
		{
			Vector<complex<T> > filterArray(N);
			for(Index i = 0; i < N; i++)
				filterArray[i] = filterData[i];
			Fourier<T>::FFT(filterArray, N, 1.);
			for(Index i = 0; i < N; i++)
				filterArray[i] *= this->fftData[i];
			Fourier<T>::IFFT(filterArray, N);
			for(int i = 0; i < N; i++)
			{
				realResult[i] = filterArray[i].real();
				if(imaginaryResult != NULL)
					imaginaryResult[i] = filterArray[i].imag();
			}
		}

		template <typename T>
			void Fourier<T>::ConvolveWith(Fourier<T>& filter, ushort N,
					Vector<T>& realResult, Vector<T>& imaginaryResult)
			{
				Vector<complex<T> > filterArray(N);
				filter.getFftData(filterArray);
				for(Index i = 0; i < N; i++)
					filterArray[i] *= N;
				for(Index i = 0; i < N; i++)
					filterArray[i] *= this->fftData[i];
				Fourier<T>::IFFT(filterArray, N);
				for(int i = 0; i < N; i++)
				{
					realResult[i] = filterArray[i].real();
					if(imaginaryResult != NULL)
						imaginaryResult[i] = filterArray[i].imag();
				}
			}

		template <typename T>
		void Fourier<T>::getData(Vector<T>& result)
		{
			for(Index i = 0; i < N; i++)
				result[i] = this->data[i];
		}

		template <typename T>
		void Fourier<T>::getFftData(Vector<T>& realPart, Vector<T>& imagPart)
		{
			for(Index i = 0; i < N; i++)
			{
				realPart[i] = this->fftData[i].real();
				if(imagPart)
					imagPart[i] = this->fftData[i].imag();
			}
		}

		template <typename T>
		void Fourier<T>::getFftData(Vector<complex<T>>& result)
		{
			for(Index i = 0; i < N; i++)
				result[i] = this->fftData[i];
		}

		template <typename T>
		void Fourier<T>::getIfftData(Vector<T>& realPart, Vector<T>& imagPart)
		{
			for(Index i = 0; i < N; i++)
			{
				realPart[i] = this->ifftData[i].real();
				if(imagPart)
					imagPart[i] = this->ifftData[i].imag();
			}
		}

		template <typename T>
		void Fourier<T>::getIfftData(Vector<complex<T> >& result)
		{
			for(Index i = 0; i < N; i++)
				result[i] = this->ifftData[i];
		}

		template <typename T>
		void Fourier<T>::FFT(Vector<T>& realPart, Vector<T>& imagPart,
				ushort N, T normalizationFactor)
		{
			Vector<complex<T>> data(N);
			for(Index i = 0; i < N; i++)
			{
				data[i].real(realPart[i]);
				if(imagPart)
					data[i].imag(imagPart[i]);
			}

			// Dapprima lo ordina col reverse order
			{
				Vector<complex<T> > temp(N);
				for(Index i = 0; i < N; i++)
					temp[i] = data[reverse(N, i)];
				for(Index j = 0; j < N; j++)
					data[j] = temp[j];
			}

			// Vettore degli zeri dell'unità.
			// Prima N/2-1 ma genera errore con ciclo for successivo in quanto prova a copiare
			// in una zona non allocata "W[N/2-1]"
			Vector<complex<T>> W(N/2);
			W[1] = polar(1., -6.28318530717959 / N);
			W[0] = 1;
			for(Index i = 2; i < N / 2; i++)
				W[i] = polar(1., -6.28318530717959 * i / N);
			int n = 1;
			int a = N / 2;
			int log2N = log2(N);
			for(Index j = 0; j < log2N; j++)
			{
				for(Index i = 0; i < N; i++)
					if(!(i & n))
					{
						// ad ogni step di raddoppiamento di n, vengono utilizzati gli indici
						// 'i' presi alternativamente a gruppetti di n, una volta si e una no.
						complex<T> temp = data[i];
						complex<T> Temp = W[(i * a) % (n * a)] * data[i + n];
						data[i] = temp + Temp;
						data[i + n] = temp - Temp;
					}
				n *= 2;
				a /= 2;
			}
			for(Index i = 0; i < N; i++)
			{
				if(normalizationFactor != 1)
					data[i] *= normalizationFactor;
				realPart[i] = data[i].real();
				if(imagPart)
					imagPart[i] = data[i].imag();
			}
		}

		template <typename T>
		void Fourier<T>::FFT(Vector<complex<T> >& data, ushort N, T normalizationFactor)
		{
			// Dapprima lo ordina col reverse order
			{
				Vector<complex<T> > temp(N);
				for(Index i = 0; i < N; i++)
					temp[i] = data[reverse(N, i)];
				for(Index j = 0; j < N; j++)
					data[j] = temp[j];
			}

			// Vettore degli zeri dell'unità.
			// Prima N/2-1 ma genera errore con ciclo for successivo in quanto prova a copiare
			// in una zona non allocata "W[N/2-1]"
			Vector<complex<T>> W(N/2);
			W[1] = polar(1., -6.28318530717959 / N);
			W[0] = 1;
			for(Index i = 2; i < N / 2; i++)
				W[i] = polar(1., -6.28318530717959 * i / N);
			int n = 1;
			int a = N / 2;
			int log2N = log2(N);
			for(Index j = 0; j < log2N; j++)
			{
				for(Index i = 0; i < N; i++)
					if(!(i & n))
					{
						// ad ogni step di raddoppiamento di n, vengono utilizzati gli indici
						// 'i' presi alternativamente a gruppetti di n, una volta si e una no.
						complex<T> temp = data[i];
						complex<T> Temp = W[(i * a) % (n * a)] * data[i + n];
						data[i] = temp + Temp;
						data[i + n] = temp - Temp;
					}
				n *= 2;
				a /= 2;
			}
			if(normalizationFactor != 1)
				for(Index i = 0; i < N; i++)
					data[i] *= normalizationFactor;
		}

		template <typename T>
		void Fourier<T>::IFFT(Vector<T>& realPart, Vector<T>& imagPart, ushort N)
		{
			Vector<complex<T> > data(N);
			for(Index i = 0; i < N; i++)
			{
				data[i].real(realPart[i]);
				if(imagPart)
					data[i].imag(imagPart[i]);
			}

			// Dapprima lo ordina col reverse order
			{
				Vector<complex<T> > temp(N);
				for(Index i = 0; i < N; i++)
					temp[i] = data[reverse(N, i)];
				for(Index j = 0; j < N; j++)
					data[j] = temp[j];
			}

			// Vettore degli zeri dell'unità.
			// Prima N/2-1 ma genera errore con ciclo for successivo in quanto prova a copiare
			// in una zona non allocata "W[N/2-1]"
			Vector<complex<T> > W(N/2);
			W[1] = polar(1., 6.28318530717959 / N);
			W[0] = 1;
			for(Index i = 2; i < N / 2; i++)
				W[i] = polar(1., 6.28318530717959 * i / N);
			int n = 1;
			int a = N / 2;
			int log2N = log2(N);
			for(Index j = 0; j < log2N; j++)
			{
				for(Index i = 0; i < N; i++)
					if(!(i & n))
					{
						// ad ogni step di raddoppiamento di n, vengono utilizzati gli indici
						// 'i' presi alternativamente a gruppetti di n, una volta si e una no.
						complex<T> temp = data[i];
						complex<T> Temp = W[(i * a) % (n * a)] * data[i + n];
						data[i] = temp + Temp;
						data[i + n] = temp - Temp;
					}
				n *= 2;
				a /= 2;
			}
			for(Index i = 0; i < N; i++)
			{
				realPart[i] = data[i].real();
				if(imagPart)
					imagPart[i] = data[i].imag();
			}
		}

		template <typename T>
		void Fourier<T>::IFFT(Vector<complex<T> >& data, ushort N)
		{
			// Dapprima lo ordina col reverse order
			{
				Vector<complex<T> > temp(N);
				for(Index i = 0; i < N; i++)
					temp[i] = data[reverse(N, i)];
				for(Index j = 0; j < N; j++)
					data[j] = temp[j];
			}

			// Vettore degli zeri dell'unità.
			// Prima N/2-1 ma genera errore con ciclo for successivo in quanto prova a copiare
			// in una zona non allocata "W[N/2-1]"
			Vector<complex<T>> W(N/2);
			W[1] = polar(1., 6.28318530717959 / N);
			W[0] = 1;
			for(Index i = 2; i < N / 2; i++)
				W[i] = polar(1., 6.28318530717959 * i / N);
			int n = 1;
			int a = N / 2;
			int log2N = log2(N);
			for(Index j = 0; j < log2N; j++)
			{
				for(Index i = 0; i < N; i++)
					if(!(i & n))
					{
						// ad ogni step di raddoppiamento di n, vengono utilizzati gli indici
						// 'i' presi alternativamente a gruppetti di n, una volta si e una no.
						complex<T> temp = data[i];
						complex<T> Temp = W[(i * a) % (n * a)] * data[i + n];
						data[i] = temp + Temp;
						data[i + n] = temp - Temp;
					}
				n *= 2;
				a /= 2;
			}
		}

		template <typename T>
		void Fourier<T>::Convolve(Vector<T>& data, Vector<T>& filter, ushort N,
				Vector<T>& realResult, Vector<T>& imaginaryResult)
		{
			Vector<complex<T>> dataArray(N);
			Vector<complex<T>> filterArray(N);
			for(Index i = 0; i < N; i++)
			{
				dataArray[i] = data[i];
				filterArray[i] = filter[i];
			}
			Fourier<T>::FFT(dataArray, N, 1. / N);
			Fourier<T>::FFT(filterArray, N, 1.);
			for(Index i = 0; i < N; i++)
				dataArray[i] *= filterArray[i];
			Fourier<T>::IFFT(dataArray, N);
			for(Index i = 0; i < N; i++)
			{
				realResult[i] = dataArray[i].real();
				if(imaginaryResult)
					imaginaryResult[i] = dataArray[i].imag();
			}
		}
	}
}
