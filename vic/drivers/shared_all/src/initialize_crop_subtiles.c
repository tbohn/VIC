/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes crop subtiles.
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
 * @brief    Initialize crop subtiles
 *****************************************************************************/
void
initialize_crop_subtiles(all_vars_struct *all_vars,
                         veg_hist_struct  *veg_hist,
                         veg_con_struct  *veg_con)
{
    extern option_struct options;

    size_t               Nveg;
    unsigned int         veg;
    unsigned int         band;
    unsigned int         i;
    double         firr;

    cell_data_struct   **cell;
    energy_bal_struct  **energy;
    snow_data_struct   **snow;
    veg_var_struct     **veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    Nveg = veg_con[0].vegetat_type_num;

    for ( veg = 0 ; veg < Nveg ; veg++ ) {
        if (options.CROPSPLIT && veg_con[veg].crop_split) {
            firr = veg_hist[veg].firr[NR];
            for (i = 0; i < veg_con[veg].Nsubtiles; i++) {
                for ( band = 0; band < options.SNOW_BAND; band++ ) {

                    copy_cell_data(&(cell[veg][band]),&(all_vars->cell_subtiles[veg][i][band]));
                    copy_veg_var(&(veg_var[veg][band]),&(all_vars->veg_var_subtiles[veg][i][band]));
                    copy_snow_data(&(snow[veg][band]),&(all_vars->snow_subtiles[veg][i][band]));
                    copy_energy_bal(&(energy[veg][band]),&(all_vars->energy_subtiles[veg][i][band]));

                    if (i == 0) {
                        all_vars->veg_var_subtiles[veg][i][band].Wdew = 0;
                    }
                    else {
                        if (firr > 0) {
                            all_vars->veg_var_subtiles[veg][i][band].Wdew = veg_var[veg][band].Wdew / firr;
                        }
                    }
                    if (i == 0) {
                        all_vars->snow_subtiles[veg][i][band].snow_canopy = 0;
                    }
                    else {
                        if (firr > 0) {
                            all_vars->snow_subtiles[veg][i][band].snow_canopy = snow[veg][band].snow_canopy / firr;
                        }
                    }

                }
            }
        }
    }
}
