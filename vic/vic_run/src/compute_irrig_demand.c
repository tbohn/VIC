/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute irrigation demand based on soil moisture.
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
 * @brief    Compute irrigation demand based on soil moisture.
 *****************************************************************************/
void
compute_irrig_demand(cell_data_struct  *cell,
                     soil_con_struct   *soil_con,
                     veg_lib_struct    *veg_lib,
                     force_data_struct *force,
                     double             swq)
{

    unsigned short thresh_idx;
    unsigned short target_idx;
    double moistfract;
    double irr_sm_thresh;
    double irr_sm_target;
    double irrig_est;

    if (swq < 0.001 && force->air_temp[NR] > 7) {

        thresh_idx = 0;
        target_idx = 0;
        moistfract = cell->layer[target_idx].moist;
        if (veg_lib->ithresh == IRR_SAT)
            irr_sm_thresh = soil_con->max_moist[thresh_idx];
        else if (veg_lib->ithresh == IRR_FC)
            irr_sm_thresh = soil_con->Wcr[thresh_idx] / 0.7;
        else
            irr_sm_thresh = soil_con->Wcr[thresh_idx]; // critical point in most cases
        if (irr_sm_thresh > soil_con->max_moist[thresh_idx])
            irr_sm_thresh = soil_con->max_moist[thresh_idx];
        if (veg_lib->itarget == IRR_SAT)
            irr_sm_target = soil_con->max_moist[target_idx];
        else
            irr_sm_target = soil_con->Wcr[target_idx] / 0.7; // field capacity in most cases

        if (irr_sm_target > soil_con->max_moist[target_idx])
            irr_sm_target = soil_con->max_moist[target_idx];
        if (moistfract < irr_sm_thresh)
            cell->irr_apply = true;
        else if (moistfract >= irr_sm_target)
            cell->irr_apply = false;

        irrig_est = soil_con->max_moist[0] - cell->layer[target_idx].moist;
        if (cell->irr_apply && force->prec[NR] < irrig_est)
          cell->irr_demand = irrig_est - force->prec[NR];
        else
          cell->irr_demand = 0;

    }

}
