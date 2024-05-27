#ifndef MyInterpolators__2D_InterpolatorBase_hpp
#define MyInterpolators__2D_InterpolatorBase_hpp

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <iterator>

namespace _2D {

    template <typename T>
    struct RealTypeOf {
        using type = double;
    };

    template <template <typename> class T, typename R>
    struct RealTypeOf<T<R>> {
        using type = R;
    };

    template <template <typename, typename> class T, typename B, typename R>
    struct RealTypeOf<T<B, R>> {
        using type = R;
    };

    template <class Derived, typename Real = typename RealTypeOf<Derived>::type>
    class InterpolatorBase {
    public:
        using VectorType = std::vector<Real>;

        std::vector<Real> xView, yView, zView;

        InterpolatorBase() = default;

        InterpolatorBase(const InterpolatorBase& rhs)
            : xView(rhs.xView), yView(rhs.yView), zView(rhs.zView) {}

        friend void swap(InterpolatorBase& lhs, InterpolatorBase& rhs) {
            using std::swap;
            swap(lhs.xView, rhs.xView);
            swap(lhs.yView, rhs.yView);
            swap(lhs.zView, rhs.zView);
        }

        InterpolatorBase& operator=(InterpolatorBase rhs) {
            swap(*this, rhs);
            return *this;
        }

        const std::vector<Real>& getXData() const { return xView; }
        const std::vector<Real>& getYData() const { return yView; }
        const std::vector<Real>& getZData() const { return zView; }
        std::vector<Real> getXData() { return xView; }
        std::vector<Real> getYData() { return yView; }
        std::vector<Real> getZData() { return zView; }

        template <typename I>
        void setUnsafeDataReference(I n, const Real* x, const Real* y, const Real* z) {
            xView.assign(x, x + n);
            yView.assign(y, y + n);
            zView.assign(z, z + n);

            callSetupInterpolator<Derived>();
        }

        template <typename XIter, typename YIter, typename ZIter>
        void setData(const XIter& x_begin, const XIter& x_end, const YIter& y_begin,
            const YIter& y_end, const ZIter& z_begin, const ZIter& z_end) {
            xView.assign(x_begin, x_end);
            yView.assign(y_begin, y_end);
            zView.assign(z_begin, z_end);

            setUnsafeDataReference(xView.size(), xView.data(), yView.data(), zView.data());
        }

        template <typename I, typename X, typename Y, typename Z>
        void setData(I n, const X* x, const Y* y, const Z* z) {
            setData(x, x + n, y, y + n, z, z + n);
        }

        template <typename I, typename X, typename Y, typename Z>
        typename std::enable_if<std::is_integral<I>::value>::type setData(
            I nx, const X* x, I ny, const Y* y, I nz, const Z* z) {
            setData(x, x + nx, y, y + ny, z, z + nz);
        }

        template <typename XT, typename YT, typename ZT>
        auto setData(const XT& x, const YT& y, const ZT& z)
            -> decltype(x.size(), x.data(), y.size(), y.data(), z.size(), z.data(), void()) {
            return setData(x.size(), x.data(), y.size(), y.data(), z.size(), z.data());
        }

        int get_x_index_to_left_of(Real x) const {
            auto it = std::lower_bound(xView.begin(), xView.end(), x);
            return std::distance(xView.begin(), it) - 1;
        }

        int get_x_index_to_right_of(Real x) const {
            return get_x_index_to_left_of(x) + 1;
        }

        int get_y_index_below(Real y) const {
            auto it = std::lower_bound(yView.begin(), yView.end(), y);
            return std::distance(yView.begin(), it) - 1;
        }

        int get_y_index_above(Real y) const {
            return get_y_index_below(y) + 1;
        }

        int get_x_index_closest_to(Real x) const {
            auto it = std::lower_bound(xView.begin(), xView.end(), x);
            auto dist_lower = std::distance(xView.begin(), it);
            if (it != xView.begin() && (it == xView.end() || std::abs(x - *(it - 1)) < std::abs(x - *it)))
                --it;
            return std::distance(xView.begin(), it);
        }

        int get_y_index_closest_to(Real y) const {
            auto it = std::lower_bound(yView.begin(), yView.end(), y);
            auto dist_lower = std::distance(yView.begin(), it);
            if (it != yView.begin() && (it == yView.end() || std::abs(y - *(it - 1)) < std::abs(y - *it)))
                --it;
            return std::distance(yView.begin(), it);
        }

    protected:
        void checkData() const {
            if (xView.empty() || yView.empty() || zView.empty())
                throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
            if (xView.size() == 0 || yView.size() == 0 || zView.size() == 0)
                throw std::logic_error("Interpolator data is zero size. Did you call setData() with non-zero sized vectors?");
        }

        template <typename T>
        struct has_setupInterpolator {
            template <typename U, void (U::*)()>
            struct SFINAE {};
            template <typename U>
            static char Test(SFINAE<U, &U::setupInterpolator>*);
            template <typename U>
            static int Test(...);
            static const bool value = sizeof(Test<T>(0)) == sizeof(char);
        };

        template <typename T>
        typename std::enable_if<has_setupInterpolator<T>::value>::type
            callSetupInterpolator() {
            static_cast<T*>(this)->setupInterpolator();
        }

        template <typename T>
        typename std::enable_if<!has_setupInterpolator<T>::value>::type
            callSetupInterpolator() {}
    };

}  // namespace _2D

#endif  // include protector
