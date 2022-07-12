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


#include <vector>
#include <iostream>
#include <cmath>

#include "AbstractEqnConc.h"

class EqnConc_3x : public AbstractEqnConc {
public:
    EqnConc_3x()
    {
        m_current_concentrations = { 0, 0, 0 };
    }

    std::vector<double> solver()
    {
        double error = 0;
        double sum_error = 0;
        double diff_error = 0;
        auto lastConcentration = m_current_concentrations;
        std::vector<double> coeffs_a(m_stoichiometry[0] + 1), coeffs_b(m_stoichiometry[1] + 1), coeffs_c(m_stoichiometry[2] + 1);


        for (m_lastIter = 0; m_lastIter < m_maxiter; ++m_lastIter) {
            std::cout << error << " " << sum_error << " " << diff_error << std::endl;
            powSpecies();

        //printConcentration(m_lastIter);

        for (int a = 1; a <= m_stoichiometry[0]; ++a) {
            coeffs_a[a] = 0;
            for (int b = 1; b <= m_stoichiometry[1]; b++) {
                for(int c = 1; c <= m_stoichiometry[2]; ++c)
                {
                    double beta = m_stability_constants[Index({ a, b, c })];
                    coeffs_a[a] += a * beta * m_powed_species[1][b - 1]* m_powed_species[1][c - 1];
                }
            }
        }

        coeffs_a[0] = -m_initial_concentrations[0];
        coeffs_a[1] += 1;
        AllConcentrations();

        lastConcentration[0] = PolynomialSolver(0, m_initial_concentrations[0], coeffs_a);

        //powSpecies();
        //PrintInfos();


        for (int b = 1; b <= m_stoichiometry[1]; b++) {
            coeffs_b[b] = 0;
            for (int a = 1; a <= m_stoichiometry[0]; ++a) {
                for(int c = 1; c <= m_stoichiometry[2]; ++c)
                {
                double beta = m_stability_constants[Index({ a, b,c })];
                coeffs_b[b] += b * beta * m_powed_species[0][a - 1]* m_powed_species[0][c - 1];
                }
            }
        }

        coeffs_b[0] = -m_initial_concentrations[1];
        coeffs_b[1] += 1;

        AllConcentrations();


        lastConcentration[1] = PolynomialSolver(0, m_initial_concentrations[1], coeffs_b);

        //powSpecies(0);
        //powSpecies(1);
       // powSpecies();

        for (int c = 1; c <= m_stoichiometry[2]; c++) {
            coeffs_c[c] = 0;
            for (int a = 1; a <= m_stoichiometry[0]; ++a) {
                for(int b = 1; b <= m_stoichiometry[1]; ++b)
                {
                double beta = m_stability_constants[Index({ a, b,c })];
                coeffs_c[c] += b * beta * m_powed_species[0][a - 1]* m_powed_species[0][b - 1];
                }
            }
        }
        coeffs_c[0] = -m_initial_concentrations[2];
        coeffs_c[1] += 1;

        lastConcentration[2] = PolynomialSolver(0, m_initial_concentrations[2], coeffs_c);


        m_lastConv = Convergency(lastConcentration);
        DampConcentration(lastConcentration, 0.5);
        //MixConcentration(lastConcentration);
        /*
        std::cout << error << std::endl;
        if(error <= 0)
            DampConcentration(lastConcentration, 0.5);
        else
            DampConcentration(lastConcentration, 0.9);
            */
        AllConcentrations();
        //Scale();
        //AllConcentrations();

        error = ConcentrationalError();
        if(diff_error < std::abs(sum_error - (sum_error + error)))
        {
            //PrintInfos();

            Scale(0.85);
            //DampConcentration(lastConcentration, 0.8);
            //Scale();
            AllConcentrations();

            std::cout << "scaled" << std::endl;
            //PrintInfos();
        }
        diff_error = std::abs(sum_error - (sum_error + error));
        sum_error += error;
        if (m_lastConv < m_converge) {
            return currentConcentration();
            }
        }
        return m_current_concentrations;
    }

    std::vector<double> AllConcentrations()
    {
        std::vector<double> vector(m_prod_species + m_current_concentrations.size());
        m_recalulated = m_current_concentrations;

        for (unsigned int i = 0; i < m_current_concentrations.size(); ++i)
            vector[i] = m_current_concentrations[i];

        int index = m_current_concentrations.size();
        for (int a = 1; a <= m_stoichiometry[0]; ++a) {
            double powA = a * pow(m_current_concentrations[0], a);
            for (int b = 1; b <= m_stoichiometry[1]; b++) {
                double powB = pow(m_current_concentrations[1], b);
                for(int c = 1; c <= m_stoichiometry[2]; c++)
                {
                    double beta = m_stability_constants[Index({ a, b, c })];
                    const double conc = (beta * powB * powA* pow(m_current_concentrations[2], c));
                    vector[index++] = conc;
                    m_recalulated[0] += a * conc;
                    m_recalulated[1] += b * conc;
                    m_recalulated[2] += c * conc;
                }
            }
        }
        return vector;
    }
};
