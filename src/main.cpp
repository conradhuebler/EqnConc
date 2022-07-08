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


#include <iostream>
#include <fstream>
#include <random>

#include "json.hpp"

#include "include/EqnConc_2.h"

using json = nlohmann::json;

int main(int argc, char **argv)
{
    EqnConc_2 solver;

    if(argc == 2)
    {
    nlohmann::json file;
       std::ifstream restart_file(argv[1]);
       std::cout << argv[1] << std::endl;
       try {
           restart_file >> file;
       } catch (nlohmann::json::type_error& e) {
           std::cout << e.what() << std::endl;
       } catch (nlohmann::json::parse_error& e) {
           std::cout << e.what() << std::endl;
       }
       std::cout << file << std::endl;

       solver.setInitialConcentrations(file["A0"].get<double>(),file["B0"].get<double>());
       solver.setStoichiometry(file["A"].get<int>(),file["B"].get<int>());
       auto v7 = file["constants"].get<std::vector<double>>();
       std::vector<double> constants;
       for (auto i : v7)
       {
           constants.push_back(pow(10, i));
           std::cout << i << " " << constants[constants.size() - 1] << std::endl;
       }

       solver.setStabilityConstants(constants);

       solver.setMaxIter(file.value("MaxIter", 1e3));
       solver.setConvergeThreshold(file.value("ConvThresh", 1e-20));
    }else
    {
        solver.setInitialConcentrations(0.00100009, 0.003979);
        solver.setStoichiometry(10,10);
        //solver.setStoichiometry(1,2);
        std::default_random_engine generator;
        std::uniform_int_distribution<int> distribution(1,6);

        std::vector<double> constants(10*10);
        for(int i = 0; i < constants.size(); ++i)
            constants[i] = distribution(generator);
        //solver.setStabilityConstants({pow(10, 3.99995), pow(10, 3.99995 + 2.49997) });
        solver.setStabilityConstants(constants);
        solver.setMaxIter(1e5);
        solver.setConvergeThreshold(1e-20);
    }
    solver.Guess();
    auto result = solver.solver();
    for(auto i : result)
        std::cout << i << " ";
    std::cout << std::endl;

    result = solver.AllConcentrations();
    for(auto i : result)
        std::cout << i << " ";
    std::cout << std::endl;

    result = solver.RecalculatedInititial();
    for(auto i : result)
        std::cout << i << " ";
    std::cout << std::endl;

    return 0;
}
