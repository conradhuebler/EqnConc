/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <cmath>
#include <iostream>
#include <vector>

class EqnConc_2 {
public:
    EqnConc_2()
    {
        m_current_concentration = { 0, 0 };
    }

    inline void setMaxIter(int maxiter)
    {
        m_maxiter = maxiter;
    }

    inline void setStabilityConstants(const std::vector<double>& vector)
    {
        m_stability_constants = vector;
    }

    inline void setConvergeThreshold(double converge)
    {
        m_converge = converge;
    }
    inline void setStoichiometry(int A, int B)
    {
        m_A = A;
        m_B = B;
    }

    inline void setInitialConcentrations(double A, double B)
    {
        m_A0 = A;
        m_B0 = B;
    }

    inline int LastIterations() const { return m_lastIter; }
    inline double LastConvergency() const { return m_lastConv; }
    std::vector<double> currentConcentration() const { return std::vector<double>() = { m_cA, m_cB }; }

    void Guess()
    {
        m_cA = std::min(m_A0, m_B0) / (10.0 * (m_A + m_B));
        m_cB = m_cA;
        m_current_concentration[0] = m_cA;
        m_current_concentration[1] = m_cB;

        m_powA.resize(m_A);
        m_powB.resize(m_B);
    }

    inline bool Converged() const { return m_converged; }
    std::vector<double> solver()
    {
        /*
        if (m_A0 < m_converge || m_B0 < m_converge)
            return std::vector<double>{ m_A0, m_B0 };
        */

        std::vector<double> coeffs_a(m_A + 1), coeffs_b(m_B + 1);
        // powA();
        powB();
        // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        // bool hasit;
        // double thresh;
        for (m_lastIter = 0; m_lastIter < m_maxiter; ++m_lastIter) {
            printConcentration(m_lastIter);
            for (int a = 1; a <= m_A; ++a) {
                coeffs_a[a] = 0;
                for (int b = 1; b <= m_B; b++) {
                    double beta = m_stability_constants[Index(a, b)];
                    // coeffs_a(a) += a * beta * pow(m_cB, b);
                    coeffs_a[a] += a * beta * m_powB[b - 1]; // pow(m_cB, b);
                }
            }
            coeffs_a[0] = -m_A0;
            coeffs_a[1] += 1;

            // coeffs_a = FillAVector();

            // if (abs(coeffs_a(coeffs_a.size() - 1) - 1) < m_converge) {
            //     solver.compute(coeffs_a);
            //     m_cA = solver.greatestRealRoot(hasit, thresh);
            // } else
            m_cA = PolynomialSolver(0, m_A0, coeffs_a);
            powA();
            /*
            m_powA.resize(0);
            m_powA.push_back(m_cA);
            for(int a = 1; a < m_A; ++a)
            {
                m_powA.push_back(m_powA[m_powA.size() -1 ]*m_cA);
            }
        */
            // coeffs_b = FillBVector();
            for (int b = 1; b <= m_B; b++) {
                coeffs_b[b] = 0;
                for (int a = 1; a <= m_A; ++a) {
                    double beta = m_stability_constants[Index(a, b)];
                    // coeffs_b(b) += b * beta * pow(m_cA, a);
                    coeffs_b[b] += b * beta * m_powA[a - 1]; // pow(m_cA, a);
                }
            }
            coeffs_b[0] = -m_B0;
            coeffs_b[1] += 1;

            // if (abs(coeffs_b(coeffs_b.size() - 1) - 1) < m_converge) {
            //     solver.compute(coeffs_b);
            //     m_cB = solver.greatestRealRoot(hasit, thresh);
            //} else
            m_cB = PolynomialSolver(0, m_B0, coeffs_b);
            powB();
            /*
                m_powB.resize(0);
                m_powB.push_back(m_cB);
                for(int a = 1; a < m_A; ++a)
                {
                    m_powB.push_back(m_powB[m_powB.size() -1 ]*m_cB);
                }
        */
            m_lastConv = std::abs(m_current_concentration[0] - m_cA) + std::abs(m_current_concentration[1] - m_cB);
            m_current_concentration[0] = m_cA;
            m_current_concentration[1] = m_cB;
            if (m_lastConv < m_converge) {
                return currentConcentration();
            }
        }
        return m_current_concentration;
    }

    inline int Index(int a, int b) const { return (a - 1) * m_B + (b - 1); }

    void printConcentration(int iter) const
    {
        std::cout << iter << " " << m_cA << " " << m_cB << std::endl;
    }
    inline void powA()
    {
        double tmp = m_cA;
        for (int a = 0; a < m_A; ++a) {
            m_powA[a] = tmp;
            tmp *= m_cA;
        }
    }

    inline void powB()
    {
        double tmp = m_cB;
        for (int b = 0; b < m_B; ++b) {
            m_powB[b] = tmp;
            tmp *= m_cB;
        }
    }
    inline int Timer() const { return m_time; }

    std::vector<double> AllConcentrations()
    {
        std::vector<double> vector(m_A * m_B + 2);
        m_recalulated = { m_current_concentration[0], m_current_concentration[1] };
        vector[0] = m_current_concentration[0];
        vector[1] = m_current_concentration[1];

        int index = 2;
        for (int a = 1; a <= m_A; ++a) {
            double powA = a * pow(m_current_concentration[0], a);
            for (int b = 1; b <= m_B; b++) {
                double beta = m_stability_constants[Index(a, b)];
                const double c = (beta * pow(m_current_concentration[1], b) * powA);
                vector[index++] = c;
                m_recalulated[0] += a * c;
                m_recalulated[1] += b * c;
            }
        }
        return vector;
    }

    std::vector<double> RecalculatedInititial() const
    {
        return m_recalulated;
    }

private:
    double PolynomialSolver(double min, double max, std::vector<double> polynom)
    {
        while (polynom[polynom.size() - 1] < m_converge)
            polynom.erase(polynom.end() - 1);

        if (polynom.size() == 2) {
            return -polynom[0] / polynom[1];
        } else if (polynom.size() == 3) {
            // if(polynom[2] < m_converge)
            //     return -polynom[0]/polynom[1];
            double p = polynom[1] / polynom[2];
            double q = polynom[0] / polynom[2];
            double p2 = p * 0.5;
            double psqaure = p2 * p2;
            double sqrt_eval = sqrt(psqaure - q);
            double x1 = -p2 + sqrt_eval;
            if (std::isnan(x1)) {
                return Bisection(min, max, polynom);
            }
            if (min <= x1 && x1 <= max)
                return x1;
            else
                return -p2 - sqrt_eval;
        } else
            return Bisection(min, max, polynom);
    }

    double Bisection(double min, double max, const std::vector<double>& polynom)
    {
        double mean = (min + max) / 2.0;
        double y_min = 0;
        double y_max = 0;
        double result = mean;

        y_min = fPolynom(polynom, min);
        y_max = fPolynom(polynom, max);

        for (int i = 0; i < m_maxiter; ++i) {
            double bisect = fPolynom(polynom, mean);
            if (std::signbit(y_min) == std::signbit(bisect)) {
                min = mean;
                y_min = fPolynom(polynom, min);
            } else if (std::signbit(y_max) == std::signbit(bisect)) {
                max = mean;
                y_max = fPolynom(polynom, max);
            }

            result = mean;
            mean = (min + max) / 2.0;
            if (std::abs(result - mean) < m_converge)
                return mean;
        }
        return result;
    }

    double fPolynom(const std::vector<double>& polynom, double x)
    {
        double result = 0.0;
        for (int i = 0; i < polynom.size(); ++i)
            result += polynom[i] * pow(x, i);
        return result;
    }

    double m_A0 = 0, m_B0 = 0;
    double m_cA = 0, m_cB = 0;
    double m_converge = 1e-10;
    int m_A = 0, m_B = 0;
    int m_maxiter = 1000;
    int m_time = 0;
    int m_lastIter = 0;
    double m_lastConv = 0.0;
    std::vector<double> m_stability_constants;
    std::vector<double> m_current_concentration;
    std::vector<double> m_recalulated;
    std::vector<double> m_powA, m_powB;
    bool m_converged = false;
};
