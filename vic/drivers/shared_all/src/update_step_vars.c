/******************************************************************************
* @section DESCRIPTION
*
* This subroutine updates data structures with values for the current
* time step.
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
******************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
* @brief        This subroutine updates data structures with values for the
*               current time step.
******************************************************************************/
int
update_step_vars(all_vars_struct *all_vars,
                 veg_con_struct  *veg_con,
                 veg_hist_struct *veg_hist)
{
    extern option_struct options;

    unsigned short       iveg;
    size_t               Nveg;
    unsigned short       band;
    size_t               Nbands;
    veg_var_struct      *veg_var;
    snow_data_struct    *snow;
    veg_var_struct      *veg_var1;
    veg_var_struct      *veg_var2;
    snow_data_struct    *snow1;
    snow_data_struct    *snow2;

    Nbands = options.SNOW_BAND;

    /* Set number of vegetation types */
    Nveg = veg_con[0].vegetat_type_num;

    /* Assign current veg characteristics */
    for (iveg = 0; iveg <= Nveg; iveg++) {
        for (band = 0; band < Nbands; band++) {

            /* set local pointers */
            veg_var = &(all_vars->veg_var[iveg][band]);
            snow = &(all_vars->snow[iveg][band]);

            veg_var->albedo = veg_hist[iveg].albedo[NR];
            veg_var->displacement = veg_hist[iveg].displacement[NR];
            veg_var->fcanopy = veg_hist[iveg].fcanopy[NR];
            veg_var->LAI = veg_hist[iveg].LAI[NR];
            veg_var->roughness = veg_hist[iveg].roughness[NR];
            if (options.IRRIGATION) {
                veg_var->fcrop = veg_hist[iveg].fcrop[NR];
                veg_var->firr = veg_hist[iveg].firr[NR];
            }

            /* Convert parameters from spatial-average to plant-specific */
            if (options.CROPSPLIT && veg_con[iveg].crop_split) {
                // Tile contains irrigated and non-irrigated sub-tiles
                veg_var1 = &(all_vars->veg_var_subtiles[iveg][0][band]);
                veg_var2 = &(all_vars->veg_var_subtiles[iveg][1][band]);
                snow1 = &(all_vars->snow_subtiles[iveg][0][band]);
                snow2 = &(all_vars->snow_subtiles[iveg][1][band]);
                convert_to_plant_specific_crop_subtiles(veg_var, veg_var1,
                                                        veg_var2, snow1,
                                                        snow2);
            }
            else {
                // Not a split tile
                convert_to_plant_specific(veg_var, snow);
            }
        }


    }

    if (options.CROPSPLIT) {
        handle_subarea_changes(all_vars, veg_con, veg_hist);
    }

    return (0);
}
