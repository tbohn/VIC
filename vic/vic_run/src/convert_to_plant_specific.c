/******************************************************************************
 * @section DESCRIPTION
 *
 * Convert to plant-specific veg parameters
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

#include <vic_run.h>

/******************************************************************************
 * @brief    Convert to plant-specific veg parameters.
 *****************************************************************************/
void
convert_to_plant_specific(veg_var_struct   *veg_var,
                          snow_data_struct *snow)
{

    extern parameters_struct param;

    if (veg_var->fcanopy > 0) {
        veg_var->LAI /= veg_var->fcanopy;
        veg_var->Wdew /= veg_var->fcanopy;
        snow->snow_canopy /= veg_var->fcanopy;
    }
    veg_var->Wdmax = veg_var->LAI * param.VEG_LAI_WATER_FACTOR;

}
