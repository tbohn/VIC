/******************************************************************************
 * @section DESCRIPTION
 *
 * Convert to plant-specific veg parameters for crop subtiles
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
 * @brief    Convert to plant-specific veg parameters for crop subtiles.
 *****************************************************************************/
void
convert_to_plant_specific_crop_subtiles(veg_var_struct   *veg_var_whole,
                                        veg_var_struct   *veg_var_sub1,
                                        veg_var_struct   *veg_var_sub2,
                                        snow_data_struct *snow_sub1,
                                        snow_data_struct *snow_sub2)
{


    extern parameters_struct param;
    double tmp_fallow_fcanopy;
    double tmp_crop_fcanopy;
    double tmp_fallow_albedo;
    double tmp_crop_albedo;
    double tmp_fallow_LAI;
    double tmp_crop_LAI;

    // Start by setting up planted and fallow parameter values
    if (veg_var_whole->fcrop > 0) {
        tmp_fallow_fcanopy = 0;
        tmp_crop_fcanopy = veg_var_whole->fcanopy / veg_var_whole->fcrop;
        tmp_fallow_albedo = veg_var_whole->albedo;
        tmp_crop_albedo = veg_var_whole->albedo;
        tmp_fallow_LAI = 0;
        tmp_crop_LAI = veg_var_whole->LAI / veg_var_whole->fcrop;
    }
    else {
        tmp_fallow_fcanopy = veg_var_whole->fcanopy;
        tmp_crop_fcanopy = veg_var_whole->fcanopy;
        tmp_fallow_albedo = veg_var_whole->albedo;
        tmp_crop_albedo = veg_var_whole->albedo;
        tmp_fallow_LAI = veg_var_whole->LAI;
        tmp_crop_LAI = veg_var_whole->LAI;
    }
    if (tmp_fallow_fcanopy > 1) tmp_fallow_fcanopy = 1;
    if (tmp_fallow_fcanopy < MIN_FCANOPY) tmp_fallow_fcanopy = 0;
    if (tmp_crop_fcanopy > 1) tmp_crop_fcanopy = 1;
    if (tmp_crop_fcanopy < MIN_FCANOPY) tmp_crop_fcanopy = 0;
    if (tmp_fallow_albedo > 1) tmp_fallow_albedo = 1;
    if (tmp_fallow_albedo < param.ALBEDO_H2O_SURF) tmp_fallow_albedo = param.ALBEDO_H2O_SURF;
    if (tmp_crop_albedo > 1) tmp_crop_albedo = 1;
    if (tmp_crop_albedo < param.ALBEDO_H2O_SURF) tmp_crop_albedo = param.ALBEDO_H2O_SURF;

    // Distribute planted and fallow parameter values across irrigated and non-irrigated sub-tiles

    // non-irrigated (planted + fallow) sub-tile
    veg_var_sub1->firr = 0;
    if (veg_var_whole->firr < 1.0) {
        veg_var_sub1->fcanopy = ( (1-veg_var_whole->fcrop)*tmp_fallow_fcanopy + (veg_var_whole->fcrop-veg_var_whole->firr)*tmp_crop_fcanopy ) / (1-veg_var_whole->firr);
        veg_var_sub1->albedo = ( (1-veg_var_whole->fcrop)*tmp_fallow_albedo + (veg_var_whole->fcrop-veg_var_whole->firr)*tmp_crop_albedo ) / (1-veg_var_whole->firr);
        veg_var_sub1->LAI = ( (1-veg_var_whole->fcrop)*tmp_fallow_LAI + (veg_var_whole->fcrop-veg_var_whole->firr)*tmp_crop_LAI ) / (1-veg_var_whole->firr);
        veg_var_sub1->Wdmax = veg_var_sub1->LAI*param.VEG_LAI_WATER_FACTOR;
        veg_var_sub1->fcrop = (veg_var_whole->fcrop - veg_var_whole->firr)/(1 - veg_var_whole->firr);
    }
    else {
        veg_var_sub1->fcanopy = 0;
        veg_var_sub1->albedo = param.ALBEDO_BARE_SOIL;
        veg_var_sub1->LAI = 0;
        veg_var_sub1->Wdmax = 0;
        veg_var_sub1->fcrop = 0;
    }
    if (veg_var_sub1->fcanopy > 0) {
        veg_var_sub1->Wdew /= veg_var_sub1->fcanopy;
        snow_sub1->snow_canopy /= veg_var_sub1->fcanopy;
    }
    else {
        veg_var_sub1->Wdew = 0;
        snow_sub1->snow_canopy = 0;
    }

    // irrigated sub-tile
    veg_var_sub2->fcanopy = tmp_crop_fcanopy;
    veg_var_sub2->albedo = tmp_crop_albedo;
    veg_var_sub2->LAI = tmp_crop_LAI;
    veg_var_sub2->Wdmax = tmp_crop_LAI*param.VEG_LAI_WATER_FACTOR;
    veg_var_sub2->firr = 1;
    veg_var_sub2->fcrop = 1;
    if (veg_var_sub2->fcanopy > 0) {
        veg_var_sub2->Wdew /= veg_var_sub2->fcanopy;
        snow_sub2->snow_canopy /= veg_var_sub2->fcanopy;
    }
    else {
        veg_var_sub2->Wdew = 0;
        snow_sub2->snow_canopy = 0;
    }

}
