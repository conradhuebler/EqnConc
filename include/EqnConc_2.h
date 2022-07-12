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

class EqnConc_2x : public AbstractEqnConc {
public:
    EqnConc_2x()
    {
        m_current_concentrations = { 0, 0 };
    }

    std::vector<double> solver()
    {
        auto lastConcentration = m_current_concentrations;
        std::vector<double> coeffs_a(m_stoichiometry[0] + 1), coeffs_b(m_stoichiometry[1] + 1);
        powSpecies(1);
        for (m_lastIter = 0; m_lastIter < m_maxiter; ++m_lastIter) {
        for (int a = 1; a <= m_stoichiometry[0]; ++a) {
            coeffs_a[a] = 0;
            for (int b = 1; b <= m_stoichiometry[1]; b++) {
                double beta = m_stability_constants[Index({ a, b })];
                coeffs_a[a] += a * beta * m_powed_species[1][b - 1];
            }
        }
        coeffs_a[0] = -m_initial_concentrations[0];
        coeffs_a[1] += 1;

        m_current_concentrations[0] = PolynomialSolver(0, m_initial_concentrations[0], coeffs_a);
        powSpecies(0);

        for (int b = 1; b <= m_stoichiometry[1]; b++) {
            coeffs_b[b] = 0;
            for (int a = 1; a <= m_stoichiometry[0]; ++a) {
                double beta = m_stability_constants[Index({ a, b })];
                coeffs_b[b] += b * beta * m_powed_species[0][a - 1];
            }
        }
        coeffs_b[0] = -m_initial_concentrations[1];
        coeffs_b[1] += 1;

        m_current_concentrations[1] = PolynomialSolver(0, m_initial_concentrations[1], coeffs_b);
        powSpecies(1);

        m_lastConv = std::abs(m_current_concentrations[0] - lastConcentration[0]) + std::abs(m_current_concentrations[1] - lastConcentration[1]);
        lastConcentration = m_current_concentrations;
        if (m_lastConv < m_converge) {
            return currentConcentration();
            }
        }
        return m_current_concentrations;
    }

    inline int Index(const std::vector<int>& indicies) const
    {
        return (indicies[0] - 1) * m_stoichiometry[1] + (indicies[1] - 1);
    }

};
