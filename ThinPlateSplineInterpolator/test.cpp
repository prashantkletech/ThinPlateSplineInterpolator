#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>
#include "ThinPlateSplineInterpolator.hpp"

using namespace Catch;
namespace _2D {
    class TestThinPlateSplineInterp : public ThinPlateSplineInterpolator<double> {};
}  // namespace _2D

void tpsOscillations()
{
    _2D::TestThinPlateSplineInterp interp;

    int nx, ny;
    double xmin, xmax, dx;
    double ymin, ymax, dy;

    nx = 10;
    ny = 5;

    xmin = -1;
    xmax = 8;

    ymin = -1;
    ymax = 3;

    dx = (xmax - xmin) / (nx - 1);
    dy = (ymax - ymin) / (ny - 1);

    _2D::ThinPlateSplineInterpolator<double>::VectorType xx(nx * ny),
        yy(nx * ny), zz(nx * ny);

    auto f = [](double x, double y) { return sin(x) * sin(y); };

    for (int i = 0; i < nx * ny; i++) {
        // gnuplot format is essentially row-major
        xx(i) = xmin + dx * (i / ny);
        yy(i) = ymin + dy * (i % ny);
        zz(i) = f(xx(i), yy(i));
    }
    interp.setData(xx, yy, zz);


    std::ofstream out;

    out.open("TPS-oscillations-input.txt");
    for (int i = 0; i < nx * ny; i++)
        out << xx(i) << " " << yy(i) << " " << zz(i) << "\n";
    out.close();

    out.open("TPS-oscillations-output.txt");
    for (int i = 0; i < 5 * nx; i++) {
        for (int j = 0; j < 5 * ny; j++) {
            out << xmin + dx * i / 5. << " " << ymin + dy * j / 5. << " "
                << interp(xmin + dx * i / 5., ymin + dy * j / 5.) << "\n";
        }
        out << "\n";
    }
    out.close();

    std::cout << (interp(0, 0) + 1 == Approx(1 + f(0, 0)).epsilon(0.00002));
    std::cout << (interp(1, 2) == Approx(f(1, 2)).epsilon(0.00002));
    std::cout << (interp(2, 1) == Approx(f(2, 1)).epsilon(0.00002));
    std::cout << (interp(2, -1) == Approx(f(2, -1)).epsilon(0.00002));
    std::cout << (interp(8, 3) == Approx(f(8, 3)).epsilon(0.00002));
    std::cout << (interp(-2, -1) == Approx(0).epsilon(0.00002));
    std::cout << (interp(10, 3) == Approx(0).epsilon(0.00002));
    std::cout << std::endl;
}


void tpsMonotonic()
{
    _2D::TestThinPlateSplineInterp interp;

    int nx, ny;
    double xmin, xmax, dx;
    double ymin, ymax, dy;

    nx = 10;
    ny = 5;

    xmin = -1;
    xmax = 8;

    ymin = -1;
    ymax = 3;

    dx = (xmax - xmin) / (nx - 1);
    dy = (ymax - ymin) / (ny - 1);

    _2D::ThinPlateSplineInterpolator<double>::VectorType xx(nx * ny),
        yy(nx * ny), zz(nx * ny);

    auto f = [](double x, double y) { return x * y + 2 * x + 3 * y; };

    for (int i = 0; i < nx * ny; i++) {
        // gnuplot format is essentially row-major
        xx(i) = xmin + dx * (i / ny);
        yy(i) = ymin + dy * (i % ny);
        zz(i) = f(xx(i), yy(i));
    }
    interp.setData(xx, yy, zz);


    {
        std::ofstream out;

        out.open("TPS-input.txt");
        for (int i = 0; i < nx * ny; i++)
            out << xx(i) << " " << yy(i) << " " << zz(i) << "\n";
        out.close();

        out.open("TPS-output.txt");
        for (int i = 0; i < 5 * nx; i++) {
            for (int j = 0; j < 5 * ny; j++) {
                out << xmin + dx * i / 5. << " " << ymin + dy * j / 5. << " "
                    << interp(xmin + dx * i / 5., ymin + dy * j / 5.) << "\n";
            }
            out << "\n";
        }
        out.close();
    }
    std::cout << (interp(0, 0) + 1 == Approx(1 + f(0, 0)).epsilon(0.00002));
    std::cout << (interp(1, 2) == Approx(f(1, 2)).epsilon(0.00002));
    std::cout << (interp(2, 1) == Approx(f(2, 1)).epsilon(0.00002));
    std::cout << (interp(2, -1) == Approx(f(2, -1)).epsilon(0.00002));
    std::cout << (interp(8, 3) == Approx(f(8, 3)).epsilon(0.00002));
    std::cout << (interp(-2, -1) == Approx(0).epsilon(0.00002));
    std::cout << (interp(10, 3) == Approx(0).epsilon(0.00002));
    std::cout << std::endl;
}

int main()
{
    tpsMonotonic();
    tpsOscillations();
    std::cout << "done";
    return 0;

}