#ifndef MyInterpolators__2D_ThinPlateSplineInterpolator_hpp
#define MyInterpolators__2D_ThinPlateSplineInterpolator_hpp

#include "MyInterpolatorBase.hpp"
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "MatrixHelper.hpp"

namespace _2D {

    template <class Real>
    class ThinPlateSplineInterpolator : public InterpolatorBase<ThinPlateSplineInterpolator<Real>, Real> {
    public:
        using BASE = InterpolatorBase<ThinPlateSplineInterpolator<Real>, Real>;
        using VectorType = std::vector<Real>;

    protected:
        using BASE::xView;
        using BASE::yView;
        using BASE::zView;

        VectorType a, b;

    public:
        template <typename I>
        ThinPlateSplineInterpolator(I n, Real* x, Real* y, Real* z) {
            this->setData(n, x, y, z);
        }

        template <typename X, typename Y, typename Z>
        ThinPlateSplineInterpolator(X& x, Y& y, Z& z) {
            this->setData(x, y, z);
        }

        ThinPlateSplineInterpolator() : BASE() {}

        ThinPlateSplineInterpolator(const ThinPlateSplineInterpolator& rhs)
            : BASE(rhs), a(rhs.a), b(rhs.b) {}

        friend void swap(ThinPlateSplineInterpolator& lhs, ThinPlateSplineInterpolator& rhs) {
            using std::swap;
            swap(lhs.a, rhs.a);
            swap(lhs.b, rhs.b);
            swap(static_cast<BASE&>(lhs), static_cast<BASE&>(rhs));
        }

        ThinPlateSplineInterpolator& operator=(ThinPlateSplineInterpolator rhs) {
            swap(*this, rhs);
            return *this;
        }

        Real operator()(Real x, Real y) const;

    protected:
        Real G(Real x, Real y, Real xi, Real yi) const;

        void setupInterpolator();
        friend BASE;
    };

    template <class Real>
    void ThinPlateSplineInterpolator<Real>::setupInterpolator() {
        size_t n = xView.size();
        a.resize(n);
        b.resize(3);

        std::vector<std::vector<Real>> M(n, std::vector<Real>(n));
        std::vector<std::vector<Real>> N(n, std::vector<Real>(3));

        for (size_t i = 0; i < n; ++i) {
            N[i][0] = 1;
            N[i][1] = xView[i];
            N[i][2] = yView[i];

            for (size_t j = 0; j < n; ++j) {
                M[i][j] = G(xView[i], yView[i], xView[j], yView[j]);
            }
        }

        std::vector<std::vector<Real>> Minv = invertMatrix(M);
        std::vector<std::vector<Real>> Ntra = transposeMatrix(N);
        std::vector<std::vector<Real>> NtraMinv = multiplyMatrices(Ntra, Minv);
        std::vector<std::vector<Real>> NtraMinvN = multiplyMatrices(NtraMinv, N);
        std::vector<std::vector<Real>> NtraMinvNinv = invertMatrix(NtraMinvN);

        VectorType z(zView.begin(), zView.end());
        VectorType NtraMinvZ = multiplyMatrixVector(NtraMinv, z);

        b = multiplyMatrixVector(NtraMinvNinv, NtraMinvZ);
        VectorType Nb = multiplyMatrixVector(N, b);

        for (size_t i = 0; i < n; ++i) {
            a[i] = 0;
            for (size_t j = 0; j < n; ++j) {
                a[i] += Minv[i][j] * (z[j] - Nb[j]);
            }
        }
    }

    template <class Real>
    Real ThinPlateSplineInterpolator<Real>::G(Real x1, Real y1, Real x2, Real y2) const {
        if (x1 == x2 && y1 == y2) return 0;

        Real r = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        return r * r * std::log(r);
    }

    template <class Real>
    Real ThinPlateSplineInterpolator<Real>::operator()(Real x, Real y) const {
        BASE::checkData();

        if (x < *std::min_element(xView.begin(), xView.end()) ||
            x > *std::max_element(xView.begin(), xView.end()) ||
            y < *std::min_element(yView.begin(), yView.end()) ||
            y > *std::max_element(yView.begin(), yView.end())) {
            return 0;
        }

        VectorType Gx(xView.size());
        for (size_t i = 0; i < xView.size(); ++i)
            Gx[i] = G(x, y, xView[i], yView[i]);

        Real f = std::inner_product(Gx.begin(), Gx.end(), a.begin(), Real(0)) + b[0] + b[1] * x + b[2] * y;

        return f;
    }

}  // namespace _2D

#endif  // include protector
