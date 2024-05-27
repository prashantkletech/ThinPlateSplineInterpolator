#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include "MyThinPlateSplineInterpolator.hpp"

namespace _2D {
    class TestThinPlateSplineInterp : public ThinPlateSplineInterpolator<double> {};
}  // namespace _2D

void tpsOscillations()
{
    _2D::TestThinPlateSplineInterp interp;

    int nx = 10, ny = 5;
    double xmin = -1, xmax = 8, dx = (xmax - xmin) / (nx - 1);
    double ymin = -1, ymax = 3, dy = (ymax - ymin) / (ny - 1);

    std::vector<double> xx(nx * ny), yy(nx * ny), zz(nx * ny);

    auto f = [](double x, double y) { return sin(x) * sin(y); };

    for (int i = 0; i < nx * ny; ++i) {
        xx[i] = xmin + dx * (i / ny);
        yy[i] = ymin + dy * (i % ny);
        zz[i] = f(xx[i], yy[i]);
    }
    interp.setData(xx, yy, zz);

    std::ofstream out;

    out.open("TPS-oscillations-input.txt");
    for (int i = 0; i < nx * ny; ++i)
        out << xx[i] << " " << yy[i] << " " << zz[i] << "\n";
    out.close();

    out.open("TPS-oscillations-output.txt");
    for (int i = 0; i < 5 * nx; ++i) {
        for (int j = 0; j < 5 * ny; ++j) {
            out << xmin + dx * i / 5. << " " << ymin + dy * j / 5. << " "
                << interp(xmin + dx * i / 5., ymin + dy * j / 5.) << "\n";
        }
        out << "\n";
    }
    out.close();
}

void tpsMonotonic()
{
    _2D::TestThinPlateSplineInterp interp;

    int nx = 10, ny = 5;
    double xmin = -1, xmax = 8, dx = (xmax - xmin) / (nx - 1);
    double ymin = -1, ymax = 3, dy = (ymax - ymin) / (ny - 1);

    std::vector<double> xx(nx * ny), yy(nx * ny), zz(nx * ny);

    auto f = [](double x, double y) { return x * y + 2 * x + 3 * y; };

    for (int i = 0; i < nx * ny; ++i) {
        xx[i] = xmin + dx * (i / ny);
        yy[i] = ymin + dy * (i % ny);
        zz[i] = f(xx[i], yy[i]);
    }
    interp.setData(xx, yy, zz);

    std::ofstream out;

    out.open("TPS-input.txt");
    for (int i = 0; i < nx * ny; ++i)
        out << xx[i] << " " << yy[i] << " " << zz[i] << "\n";
    out.close();

    out.open("TPS-output.txt");
    for (int i = 0; i < 5 * nx; ++i) {
        for (int j = 0; j < 5 * ny; ++j) {
            out << xmin + dx * i / 5. << " " << ymin + dy * j / 5. << " "
                << interp(xmin + dx * i / 5., ymin + dy * j / 5.) << "\n";
        }
        out << "\n";
    }
    out.close();
}

int main()
{
    tpsMonotonic();
    tpsOscillations();
    std::cout << "done" << std::endl;
    return 0;
}
