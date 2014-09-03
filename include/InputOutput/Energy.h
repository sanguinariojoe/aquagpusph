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

#ifndef REPORT_ENERGY_H_INCLUDED
#define REPORT_ENERGY_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Report.h>

namespace Aqua{
namespace InputOutput{

/** \class Energy Energy.h InputOutput/Energy.h
 * Energy file loader/saver. The energy file is just a tabulated ASCII file
 * where the energy components are printed by rows (a time instant per line).
 * The energy components printed are:
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
 */
class Energy : public Report
{
public:
    /** Constructor
     */
    Energy();

    /** Destructor
     */
    ~Energy();

    /** Save the data. The data
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** Useless method
     * @return false.
     */
    bool load(){return false;}

    /** Get the file handler
     * @return The output file handler.
     */
    FILE* fileHandler(){return _file;}
private:
    /** Create the output file
     * @return false if all gone right, true otherwise.
     */
    bool create();

    /** Close the output file
     * @return false if all gone right, true otherwise.
     */
    bool close();

    /// Output file
    FILE *_file;

};  // class InputOutput

}}  // namespaces

#endif // REPORT_ENERGY_H_INCLUDED
