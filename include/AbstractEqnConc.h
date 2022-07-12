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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

class AbstractEqnConc {
public:
    AbstractEqnConc() = default;

    ~AbstractEqnConc() = default;

    inline void setMaxIter(int maxiter)
    {
        m_maxiter = maxiter;
    }

    inline void setStabilityConstants(const std::vector<double>& stability_constants)
    {
        m_stability_constants = stability_constants;
    }

    inline void setConvergeThreshold(double converge)
    {
        m_converge = converge;
    }
    inline void setStoichiometry(const std::vector<int>& stoichiometry)
    {
        m_stoichiometry = stoichiometry;
        m_prod_species = 1;
        m_sum_species = 0;
        int i = 0;
        for (int species : m_stoichiometry) {
            m_prod_species *= species;
            m_sum_species += species;
            std::vector<double> powed(m_stoichiometry[i]);
            m_powed_species.push_back(powed);
            ++i;
        }
    }

    inline void setInitialConcentrations(const std::vector<double>& initial_concentrations)
    {
        m_initial_concentrations = initial_concentrations;
    }

    inline int LastIterations() const
    {
        return m_lastIter;
    }

    inline double LastConvergency() const
    {
        return m_lastConv;
    }

    std::vector<double> currentConcentration() const
    {
        return m_current_concentrations;
    }

    virtual void Guess()
    {
        std::vector< int >blacklist;
        double min = m_initial_concentrations[0];
        for(unsigned int i = 1; i < m_initial_concentrations.size(); ++i)
        {
            if(m_initial_concentrations[i] < m_converge)
            {
                blacklist.push_back(i);
                continue;
            }
            min = std::min(min, m_initial_concentrations[i]);
        }

        min /= (10.0 * (m_sum_species));
        for(unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            m_current_concentrations[i] = min;
        for(int i : blacklist)
            m_current_concentrations[i] = 0;
    }

    inline bool Converged() const
    {
        return m_converged;
    }

    virtual std::vector<double> solver() = 0;

    virtual int Index(const std::vector<int>& indicies) const
    {
        int index = 0;
        for (unsigned int i = 0; i < indicies.size() - 1; ++i)
            index += (indicies[i] - 1) * m_stoichiometry[i + 1];
        index += indicies[indicies.size() - 1] - 1;
        return index;
    }

    void printConcentration(int iter) const
    {
        std::cout << iter;
        for (double c : m_current_concentrations)
            std::cout << " " << c;
        std::cout << std::endl;
    }

    void powSpecies(int species)
    {
        double tmp = m_current_concentrations[species];
        for (int index = 0; index < m_stoichiometry[species]; ++index) {
            m_powed_species[species][index] = tmp;
            tmp *= m_current_concentrations[species];
        }
    }

    inline int Timer() const { return m_time; }

    virtual std::vector<double> AllConcentrations()
    {
        std::vector<double> vector(m_prod_species + m_current_concentrations.size());
        m_recalulated = { m_current_concentrations[0], m_current_concentrations[1] };

        for (unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            vector[i] = m_current_concentrations[i];

        int index = m_current_concentrations.size();
        for (int a = 1; a <= m_stoichiometry[0]; ++a) {
            double powA = a * pow(m_current_concentrations[0], a);
            for (int b = 1; b <= m_stoichiometry[1]; b++) {
                double beta = m_stability_constants[Index({ a, b })];
                const double c = (beta * pow(m_current_concentrations[1], b) * powA);
                vector[index++] = c;
                m_recalulated[0] += a * c;
                m_recalulated[1] += b * c;
            }
        }
        return vector;
    }

    void UpdateConcentration(const std::vector<double> &current)
    {
        for(unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            m_current_concentrations[i] = current[i];
    }

    void MixConcentration(const std::vector<double> &current)
    {
        for(unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            m_current_concentrations[i] = (m_current_concentrations[i] + current[i])/2.0;
    }

    void DampConcentration(const std::vector<double> &current, double value)
    {
        for(unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            m_current_concentrations[i] = ((1-value)*m_current_concentrations[i] + value*current[i]);
    }

    void Scale(double value = 1)
    {
        for(unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            m_current_concentrations[i] /= m_recalulated[i]/m_initial_concentrations[i]*value;
    }

    double ConcentrationalError() const
    {
        double error = 0;
        for(int i = 0; i < m_initial_concentrations.size(); ++i)
        {
            error +=( m_initial_concentrations[i] - m_recalulated[i]);
        }
        return error;
    }

    std::vector<double> RecalculatedInititial() const
    {
        return m_recalulated;
    }

    double Convergency(const std::vector<double> &concentrations)
    {
        double result = 0.0;
        for(unsigned int i = 0; i < concentrations.size(); ++i)
        {
            result += std::abs(concentrations[i] - m_current_concentrations[i]);
        }
        return result;
    }
    void powSpecies()
    {
        for(unsigned int i = 0; i < m_initial_concentrations.size(); ++i)
            powSpecies(i);
    }

    void PrintInfos()
    {
        std::cout << "Result of Calculation - equilibrium concentration of pure species" << std::endl;
        for(auto i : m_current_concentrations)
            std::cout << i << " ";
        std::cout << std::endl << std::endl;

        std::cout << "Result of Calculation - equilibrium concentration of all species" << std::endl;
        auto result = AllConcentrations();
        for(auto i : result)
            std::cout << i << " ";
        std::cout << std::endl << std::endl;

        std::cout << "Recalculated initial concentrations" << std::endl;
        result = RecalculatedInititial();
        for(auto i : result)
            std::cout << i << " ";
        std::cout << std::endl;
    }
private:


    double NewtonRoot(double min, double max, const std::vector<double>& polynom)
    {
        double x = (min +  max)/2.0;
        std::pair<double, double> y(0,0);
        for(int iter = 0; iter < m_maxiter; ++iter)
        {
           y = function(polynom, x);
           if(std::abs(y.first) < m_converge)
               break;
           x = x - y.first/y.second;
        }
        return x;
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

    std::pair<double, double> function(const std::vector<double>& polynom, double x)
    {
        std::pair<double, double> result(0,0);
        result.first = polynom[0];
        for (unsigned int i = 1; i < polynom.size(); ++i)
        {
            result.first += polynom[i] * pow(x, i);
            result.second += i*polynom[i] * pow(x, i-1);
        }
        return result;
    }
    double fPolynom(const std::vector<double>& polynom, double x)
    {
        double result = 0.0;
        for (unsigned int i = 0; i < polynom.size(); ++i)
            result += polynom[i] * pow(x, i);
        return result;
    }

protected:
    double PolynomialSolver(double min, double max, std::vector<double> polynom)
    {
        if(std::abs(min - max) < m_converge)
            return min;

        while (std::abs(polynom[polynom.size() - 1]) < m_converge || std::isnan(polynom[polynom.size() - 1]))
            polynom.erase(polynom.end() - 1);

        if (polynom.size() < 2)
            return 0;

        if (polynom.size() == 2) {
            return -polynom[0] / polynom[1];
        } else if (polynom.size() == 3) {
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
        {
            double val = NewtonRoot(min, max, polynom);
            if(std::isnan(val))
            {
                return Bisection(min, max, polynom);
            }
            return val;
        }
    }

    std::vector<std::vector<double>> m_powed_species;

    std::vector<double> m_stability_constants;
    std::vector<double> m_current_concentrations;
    std::vector<double> m_initial_concentrations;
    std::vector<double> m_recalulated;
    std::vector<int> m_stoichiometry;

    double m_converge = 1e-10;
    double m_lastConv = 0.0;

    int m_prod_species = 0, m_sum_species;
    int m_maxiter = 1000;
    int m_time = 0;
    int m_lastIter = 0;

    bool m_converged = false;
};
