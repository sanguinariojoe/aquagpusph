/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @brief Fluid energy report.
 * (See Aqua::InputOutput::Energy for details)
 */

#ifndef REPORT_ENERGY_H_INCLUDED
#define REPORT_ENERGY_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Report.h>

namespace Aqua{
namespace InputOutput{

/** @class Energy Energy.h InputOutput/Energy.h
 * @brief Fluid energy report saver.
 *
 * In this report the fluid energy components are saved.
 *
 * This report is a plain text file.
 * Each line corresponds to a different time instant.
 * At each line the following fields are saved, separated by tabulators:
 *   -# Potential energy: \f$ E_{pot} = - \sum_i m_i
     \mathbf{g} \cdot \mathbf{r}_i \f$.
 *   -# Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
     \vert \mathbf{u}_i \vert^2 \f$.
 *   -# Internal energy: \f$ U = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \left(
        \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
        - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
     \right) m_i \mathrm{d}t \f$.
 *   -# Enthalpy: \f$ H = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
 *   -# Entropy: \f$ TS = U - H \f$.
 *   -# Total energy: \f$ E = U + E_{kin} \f$.
 *
 * The output file will be the first non existent file called `"energy.%d.dat"`,
 * where `"%d"` is replaced by a unsigned integer.
 */
class Energy : public Report
{
public:
    /// Constructor
    Energy();

    /// Destructor
    ~Energy();

    /** @brief Save the data.
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** @brief Get the energy report file handler.
     * @return The output file handler.
     */
    FILE* fileHandler(){return _file;}
private:
    /** @brief Create the output file
     * @return false if all gone right, true otherwise.
     */
    bool create();

    /** @brief Close the output file
     * @return false if all gone right, true otherwise.
     */
    bool close();

    /// Output file
    FILE *_file;

};  // class InputOutput

}}  // namespaces

#endif // REPORT_ENERGY_H_INCLUDED
