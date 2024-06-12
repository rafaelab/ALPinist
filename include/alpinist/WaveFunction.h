#ifndef ALPINIST_WAVEFUNCTION_H
#define ALPINIST_WAVEFUNCTION_H

#include <cmath>
#include <complex>
#include <type_traits>

#include <crpropa/Common.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "alpinist/CandidateProperties.h"
#include "alpinist/PlasmaDensity.h"
#include "alpinist/Constants.h"



namespace alpinist {


template<size_t N, typename T>
inline Eigen::Matrix<T, N, N> getCanonicalBases() {
	return Eigen::Matrix<T, N, N>::Identity();
}


/**
 @class WaveFunction
 @brief Wave function describing a mixed state.
  It can be written as:
   |\Psi> = sum_k c_k |\psi_k>
  Bases can be complex or real, according to the template. 
  Coefficients are always complex.
 @todo - Optimise return types to avoid needless treatment as complex.
 */
template <size_t N, class T>
class WaveFunction {
	protected:
		Eigen::Matrix<T, N, N> bases; // columns are each basis element
		Eigen::RowVector<T, N> coefficients;
	public:
		WaveFunction() {
			setBases(getCanonicalBases<N, T>());
		}

		WaveFunction(Eigen::RowVector<T, N> coeff) {
			setCoefficients(coeff);
			setBases(getCanonicalBases<N, T>());
		}

		WaveFunction(Eigen::Vector<T, N> coeff) {
			setCoefficients(coeff.transpose());
			setBases(getCanonicalBases<N, T>());
		}

		WaveFunction(Eigen::RowVector<T, N> coeff, Eigen::Matrix<T, N, N> basisSet) {
			setCoefficients(coeff);
			setBases(basisSet);
		}

		WaveFunction(Eigen::Vector<T, N> coeff, Eigen::Matrix<T, N, N> basisSet) {
			coeff.transposeInPlace();
			setCoefficients(coeff);
			setBases(basisSet);
		}

		void setBases(Eigen::Matrix<T, N, N> basisSet) {
			bases = basisSet;
		}

		void setCoefficients(Eigen::RowVector<T, N> coeff) {
			coefficients = coeff;
			normaliseCoefficients();
		}

		void normaliseCoefficients() {
			double norm = 0;
			for (size_t i = 0; i < N; i++) {
				norm += std::norm(coefficients(i));
			}
			coefficients /= norm;
		}

		void setBasis(const size_t& index, Eigen::Vector<T, N> basis) {
			bases(index, 0) = basis(0);
			bases(index, 1) = basis(1);
			bases(index, 2) = basis(2);
		}

		void setCoefficient(const size_t& index, T coeff) {
			coefficients(index) = coeff;
		}

		Eigen::Matrix<T, N, N> getBases() const {
			return bases;
		}

		Eigen::RowVector<T, N> getCoefficients() const {
			return coefficients;
		}

		Eigen::Vector<T, N> getBasis(size_t index) const {
			return bases.col(index);
		}

		T getCoefficient(size_t index) const {
			return coefficients(index);
		}

		Eigen::Vector<T, N> getKet() const {
			Eigen::Vector<T, N> v;
			for (size_t i = 0; i < N; i++) {
				v += getCoefficient(i) * getBasis(i);
			}
			return v;
		}

		Eigen::RowVector<T, N> getBra() const {
			Eigen::Vector<T, N> v;
			for (size_t i = 0; i < N; i++) {
				v += std::conj(getCoefficient(i)) * getBasis(i);
			}
			return v;
		}
		
		Eigen::RowVector<T, N> getDual() const {
			return getBra();
		}

		bool haveSameBases(const WaveFunction<N, T>& psi) const {
			bool status = false;
			for (size_t i = 0; i < N; i++) {
				if (getBasis(i) == psi.getBasis(i))
					status = true;
			}
			return status;
		}

		T innerProduct(const WaveFunction<N, T>& psi) const {
			if (haveSameBases(psi)) {
				return coefficients.conjugate().dot(psi.getCoefficients());
			}
			return getBra().dot(psi.getKet());
		}

		Eigen::Matrix<T, N, N> outerProduct(const WaveFunction<N, T>& psi) const {
			if (haveSameBases(psi)) {
				return coefficients.transpose() * psi.getCoefficients().conjugate();
			}
			return getKet() * psi.getBra();
		}

		void normalise() {
			normaliseCoefficients();
		}

		Eigen::Matrix<T, N, N> getDensityMatrix() const {
			return outerProduct(*this);
		}

		WaveFunction<N, T> selectComponent(const size_t& index) const {
			// NOTE: not a proper wavefunction because it is not normalised!
			WaveFunction<N, T> psi = *this;
			for (size_t i = 0; i < index; i++) {
				if (i != index) {
					psi.setCoefficient(i, 0.);
				}
			}

			return psi;
		}

		void operate(const Eigen::Matrix<T, N, N>& O) {
			coefficients = (O * coefficients.transpose()).transpose();
		}

		void operator+=(const WaveFunction<N, T>& other) {
			if (haveSameBases(other)) {
				coefficients += other.getCoefficients();
			} else {
				throw std::runtime_error("Cannot operate two wave functions with different bases (yet).");
			}
			normaliseCoefficients();
		}

		void operator-=(const WaveFunction<N, T>& other) const {
			if (haveSameBases(other)) {
				coefficients -= other.getCoefficients();
			} else {
				throw std::runtime_error("Cannot operate two wave functions with different bases (yet).");
			}
			normaliseCoefficients();
		}

		WaveFunction<N, T> operator+(const WaveFunction<N, T>& other) const {
			WaveFunction<N, T> psi = *this;
			psi += other;
			return psi;
		}

		WaveFunction<N, T> operator-(const WaveFunction<N, T>& other) const {
			WaveFunction<N, T> psi = *this;
			psi -= other;
			return psi;
		}

		T operator[](size_t index) const {
			return getCoefficient(index);
		}

		bool operator==(const WaveFunction<N, T>& other) const {
			if (not haveSameBases(other)) 
				return false;

			if (coefficients == other.getCoefficients())
				return true;

			return false;
		}
};


// Convenient aliases
typedef WaveFunction<3, std::complex<double>> WaveFunction3c;
typedef WaveFunction<2, std::complex<double>> WaveFunction2c;
typedef WaveFunction<1, std::complex<double>> WaveFunction1c;
typedef WaveFunction<3, double> WaveFunction3d;
typedef WaveFunction<2, double> WaveFunction2d;
typedef WaveFunction<1, double> WaveFunction1d;






} // namespace alpinist

#endif

