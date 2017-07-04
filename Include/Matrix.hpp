/*
 * Matrix.hpp
 *
 *  Created on: 20 giu 2017
 *      Author: enzo
 */

#pragma once

#include <vector>
#include "Maintainance.hpp"
#include <cassert>

using namespace std;

namespace Library {
    namespace Data {
        template<typename T>
        class Matrix {
            public:
                Matrix() {
                }
                Matrix(ushort_t rows, ushort_t cols);
                Matrix(Matrix<T>* other);
                Matrix(Matrix<T>& other);

                vector<T>& operator [](ushort_t i) {
                    return data[i];
                }
                vector<T> operator *(vector<T>& other);

                Matrix<T> inverse();
                ushort_t getRows() {
                    return data.size();
                }
                ushort_t getCols() {
                    return data[0].size();
                }
                vector<T> solveLinear(vector<T>& b);

                void fromMatrix(Matrix<T> matrix);
            private:
                T Divide(T a, T b);
                T Sign(T a);
                vector<vector<T> > data;
        };

        template<typename T>
        Matrix<T>::Matrix(ushort_t rows, ushort_t cols) :
                data(rows) {
            for (ushort_t i = 0; i < rows; i++)
                data[i] = vector<T>(cols);
        }

        template<typename T>
        Matrix<T>::Matrix(Matrix<T>& other) :
                data(other.getRows()) {
            for (ushort_t i = 0; i < other.getRows(); i++) {
                data[i] = vector<T>(other.getCols());
                for (ushort_t j = 0; j < other.getCols(); j++)
                    data[i][j] = other[i][j];
            }
        }

        template<typename T>
        Matrix<T>::Matrix(Matrix<T>* other) :
                data(other->getRows()) {
            for (ushort_t i = 0; i < other->getRows(); i++) {
                data[i] = vector<T>(other->getCols());
                for (ushort_t j = 0; j < other->getCols(); j++)
                    data[i][j] = (*other)[i][j];
            }
        }

        template<typename T>
        vector<T> Matrix<T>::operator*(vector<T>& other) {
            vector<T> result(getRows());
            for (ushort_t i = 0; i < getRows(); i++) {
                result[i] = 0;
                for (ushort_t j = 0; j < getCols(); j++)
                    result[i] += data[i][j] * other[j];
            }
            return result;
        }

        template<typename T>
        Matrix<T> Matrix<T>::inverse() {
            Matrix<T> result(getRows(), 2 * getCols());
            ushort_t rank = getRows();
            // *** Inizializzazione della matrice B = [A|I] *** //
            for (ushort_t i = 0; i < rank; i++) {
                for (ushort_t j = 0; j < rank; j++)	// *** Popolo la meta' sinistra di B con la matrice A *** //
                    result[i][j] = data[i][j];
                for (ushort_t j = rank; j < 2 * rank; j++)
                    result[i][j] = (int) (j - rank == i);// *** Popolo la meta' destra di B con la matrice identita' *** //
            }
            for (ushort_t k = 0; k < rank; k++)	// *** Azzero tutti gli elementi fuori diagonale *** //
                    {
                if (result[k][k] == 0)// *** L'elemento sulla diagonale non pu� mai essere zero! *** //
                    for (ushort_t i = k + 1; i < rank; i++)
                        if (result[i][k] != 0)// *** Scambia la i-esima riga con la k-esima *** //
                                {
                            for (ushort_t j = 0; j < 2 * rank; j++)
                                std::swap(result[i][j], result[k][j]);
                            break;// *** Fatto. Passa al calcolo dell'inversa *** //
                        }

                if (result[k][k] != 1)// *** Divido in modo da avere l'elemento sulla diagonale pari a 1*** //
                    for (ushort_t j = 2 * rank; j > k; j--)
                        result[k][j - 1] = Divide(result[k][j - 1],
                                result[k][k]);// *** Divisione stabilizzata *** //
                for (ushort_t i = 0; i < rank; i++)	// *** Sottraggo la k-esima riga alla i-esima *** //
                    if (k != i)
                        for (ushort_t j = 2 * rank; j > k; j--)	// *** Restringo agli elementi che stanno sulla diagonale *** //
                            result[i][j - 1] -= result[i][k] * result[k][j - 1];// *** Quelli a sx di essa sono gi� nulli a questo punto *** //
            }
            for (ushort_t i = 0; i < rank; i++)	// *** Copio l'inversa in A *** //
                for (ushort_t j = rank; j < 2 * rank; j++)
                    data[i][j - rank] = result[i][j];
            return *this;
        }

        template<typename T>
        T Matrix<T>::Divide(T a, T b) {
            if (a == 0 && b != 0)
                return 0;
            if (a != 0 && b == 0)
                return Sign(a) * 1E6;
            if (a == 0 && b == 0)
                return 1; // The Scream
            if (b == 1)
                return a;
            if (b == -1)
                return -1 * a;
//	if(std::abs(b) < 1E-10) return Sign(b) * a * 1E+10;

//	return Sign(a) * Sign(b) * exp(log(std::abs(a)) - log(std::abs(b)));
            return a / b;
        }

        template<typename T>
        T Matrix<T>::Sign(T a) {
            if (a < 0)
                return -1;
            return 1;
        }

        template<typename T>
        void Matrix<T>::fromMatrix(Matrix<T> matrix) {
            for (ushort_t row = 0; row < getRows(); row++)
                for (ushort_t col = 0; col < getCols(); col++)
                    data[row][col] = matrix[row][col];
        }

        template<typename T>
        vector<T> Matrix<T>::solveLinear(vector<T>& b) {
            Matrix<T> result(getRows(), getCols() + 1);
            ushort_t rank = getRows();
            vector<ushort_t> solutionIndexes(rank);
            for (ushort_t row = 0; row < rank; row++) {
                solutionIndexes[row] = row;
                result[row][rank] = b[row];
                for (ushort_t col = 0; col < rank; col++)
                    result[row][col] = (*this)[row][col];
            }
            for (ushort_t row = 0; row < rank; row++) {
                // remove zeros from the diagonal
                if (result[row][row] == 0) {
                    for (ushort_t next = row + 1; next < rank; next++) {
                        if (result[next][row] != 0) {
                            for (ushort_t col = 0; col <= rank; col++)
                                std::swap(result[next][col], result[row][col]);
                            std::swap(solutionIndexes[row],
                                    solutionIndexes[next]);
                        }
                    }
                }
                //pivoting
                if (result[row][row] != 1)
                    for (short col = rank; col >= 0; col--)
                        result[row][col] = Divide(result[row][col],
                                result[row][row]);

                //let's solve this madafakkah
                for (ushort_t otherRow = row + 1; otherRow < rank; otherRow++) {
                    if (otherRow != row) {
                        T factor = result[otherRow][row];
                        for (short col = rank; col >= row; col--) {
                            T element = Divide(result[otherRow][col], factor);
                            result[otherRow][col] = element - result[row][col];
                        }
//                        result[otherRow][row] = 0;
                    }
                }
            }
//            for (short row = rank - 1; row >= 0; row--) {
//                for (short otherRow = row - 1; otherRow >= 0; otherRow--) {
//                    T factor = (T) result[otherRow][row];
//                    for (short col = rank; col >= row; col--) {
//                        T element = Divide(result[otherRow][col], factor);
//                        result[otherRow][col] = element - result[row][col];
//                    }
//                }
//            }
            // extract the solution
            vector<T> solution(rank);
            for(short row = rank - 1; row >= 0; row--)
            {
                solution[solutionIndexes[row]] = result[row][rank];
                for(ushort_t col = row + 1; col < rank; col++)
                    solution[solutionIndexes[row]] -= result[row][col] *
                        solution[solutionIndexes[col]];
            }
            return solution;
        }
    }
}
