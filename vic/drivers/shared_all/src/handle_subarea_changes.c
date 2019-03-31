/******************************************************************************
 * @section DESCRIPTION
 *
 * Account for changes in sub-tile areas
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
 * @brief    Account for changes in sub-tile areas.
 *****************************************************************************/
void
handle_subarea_changes(all_vars_struct *all_vars,
                       veg_con_struct  *veg_con,
                       veg_hist_struct *veg_hist)
{
    extern option_struct options;

    unsigned short       iveg;
    size_t               Nveg;
    unsigned short       band;
    size_t               Nbands;
    unsigned short       l;
    veg_var_struct      *veg_var;
    cell_data_struct    *cell1;
    veg_var_struct      *veg_var1;
    snow_data_struct    *snow1;
    cell_data_struct    *cell2;
    veg_var_struct      *veg_var2;
    snow_data_struct    *snow2;
    double               old_firr;
    double               new_firr;

    Nbands = options.SNOW_BAND;

    /* Set number of vegetation types */
    Nveg = veg_con[0].vegetat_type_num;

    for(iveg = 0; iveg <= Nveg; iveg++){

        if (options.CROPSPLIT && veg_con[iveg].crop_split && veg_con[iveg].Cv > 0.0) {

            Nbands = options.SNOW_BAND;
            for ( band = 0; band < Nbands; band++ ) {

                /* set local pointers */
                veg_var = &(all_vars->veg_var[iveg][band]);
                cell1 = &(all_vars->cell_subtiles[iveg][0][band]);
                veg_var1 = &(all_vars->veg_var_subtiles[iveg][0][band]);
                snow1 = &(all_vars->snow_subtiles[iveg][0][band]);
                cell2 = &(all_vars->cell_subtiles[iveg][1][band]);
                veg_var2 = &(all_vars->veg_var_subtiles[iveg][1][band]);
                snow2 = &(all_vars->snow_subtiles[iveg][1][band]);
                if (veg_var->firr_save < 0)
                    old_firr = veg_hist[iveg].firr[0];
                else
                    old_firr = veg_var->firr_save;
                new_firr = veg_hist[iveg].firr[0];
                if (new_firr != old_firr) {
                    // The portion that grows needs to absorb state
                    // variables from the portion that shrinks; canopy
                    // storages need to be rescaled
                    if (new_firr > old_firr) { // crop fraction is growing
                        for(l = 0; l < options.Nlayer ; l++) {
                            cell2->layer[l].moist = ( cell2->layer[l].moist * old_firr + cell1->layer[l].moist * (new_firr - old_firr) ) / new_firr;
                        }
                        veg_var2->Wdew = ( veg_var2->Wdew * old_firr + veg_var1->Wdew * (new_firr - old_firr) ) / new_firr;
                        snow2->swq = ( snow2->swq * old_firr + snow1->swq * (new_firr - old_firr) ) / new_firr;
                        snow2->snow_canopy = ( snow2->snow_canopy * old_firr + snow1->snow_canopy * (new_firr - old_firr) ) / new_firr;
                    }
                    else { // fallow fraction is growing
                        for(l=0;l<options.Nlayer;l++) {
                            cell1->layer[l].moist = ( cell1->layer[l].moist * (1 - old_firr) + cell2->layer[l].moist * (old_firr - new_firr) ) / (1 - new_firr);
                        }
                        veg_var1->Wdew = ( veg_var1->Wdew * (1 - old_firr) + veg_var2->Wdew * (old_firr - new_firr) ) / (1 - new_firr);
                        snow1->swq = ( snow1->swq * (1 - old_firr) + snow2->swq * (old_firr - new_firr) ) / (1 - new_firr);
                        snow1->snow_canopy = ( snow1->snow_canopy * (1 - old_firr) + snow2->snow_canopy * (old_firr - new_firr) ) / (1 - new_firr);
                    }
                }
                veg_var->firr_save = new_firr;
            }
        }
    }
}
