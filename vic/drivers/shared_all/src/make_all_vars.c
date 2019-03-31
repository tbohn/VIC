/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine creates an array of structures that contain information about a
 * cell's states and fluxes.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Creates an array of structures that contain information about a
 *           cell's states and fluxes.
 *****************************************************************************/
all_vars_struct
make_all_vars(size_t nveg)
{
    extern option_struct     options;

    all_vars_struct temp;
    size_t          Nitems;
    size_t          v;

    Nitems = nveg + 1;

    temp.snow = make_snow_data(Nitems);
    temp.energy = make_energy_bal(Nitems);
    temp.veg_var = make_veg_var(Nitems);
    temp.cell = make_cell_data(Nitems);

    if (options.CROPSPLIT) {

        temp.snow_subtiles = (snow_data_struct ***) calloc(Nitems, sizeof(snow_data_struct **));
        temp.energy_subtiles = (energy_bal_struct ***) calloc(Nitems, sizeof(energy_bal_struct **));
        temp.veg_var_subtiles = (veg_var_struct ***) calloc(Nitems, sizeof(veg_var_struct **));
        temp.cell_subtiles = (cell_data_struct ***) calloc(Nitems, sizeof(cell_data_struct **));

        for (v = 0; v < Nitems; v++) {
            temp.snow_subtiles[v]   = make_snow_data(MAX_SUBTILES);
            temp.energy_subtiles[v] = make_energy_bal(MAX_SUBTILES);
            temp.veg_var_subtiles[v]  = make_veg_var(MAX_SUBTILES);
            temp.cell_subtiles[v]     = make_cell_data(MAX_SUBTILES);
        }

    }

    return (temp);
}
