/******************************************************************************
 * @section DESCRIPTION
 *
 * ARNO Model of Evaporation
 *
 * Routine to compute evaporation based on the assumption that evaporation is
 * at the potential for the area which is saturated, and at some percentage of
 * the potential for the area which is partial saturated.
 *
 * Evaporation from bare soil calculated only from uppermost layer.
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
 * @brief    Calculate Evaporation using ARNO model.
 *****************************************************************************/
double
arno_evap(layer_data_struct *layer,
          double             rad,
          double             air_temp,
          double             vpd,
          double            *max_moist,
          double             elevation,
          double             b_infilt,
          double             ra,
          double             delta_t,
          double            *resid_moist,
          double            *frost_fract)
{
    extern parameters_struct param;
    extern option_struct     options;

    int                      num_term;
    int                      i;
    size_t                   frost_area;
    size_t                   lidx;
    size_t                   Nlayers_esoil;
    double                   tmp, beta_asp, dummy;
    double                   ratio, as;
    double                   Epot; /* potential bare soil evaporation */
    double                   moist;
    double                   esoil;
    double                   max_infil;
    double                   Evap;
    double                   tmpsum;
    double                   tmpesoil;
    double                   max_moist_tot;
    double                   resid_moist_tot;

    Evap = 0;

    /* moist = liquid soil moisture */
    moist = 0;
    max_moist_tot = 0;
    resid_moist_tot = 0;
    if (options.DEEP_ESOIL) {
        Nlayers_esoil = 2;
    }
    else {
        Nlayers_esoil = 1;
    }
    for (lidx = 0; lidx < Nlayers_esoil; lidx++) {
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            moist +=
                (layer[lidx].moist -
                 layer[lidx].ice[frost_area]) * frost_fract[frost_area];
        }
        max_moist_tot += max_moist[lidx];
        resid_moist_tot += resid_moist[lidx];
    }
    if (moist > max_moist_tot) {
        moist = max_moist_tot;
    }

    /* Calculate the potential bare soil evaporation (mm/time step) */

    Epot =
        penman(air_temp, elevation, rad, vpd, ra, 0.0,
               param.SOIL_RARC) * delta_t / CONST_CDAY;

    /**********************************************************************/
    /*  Compute temporary infiltration rate based on given soil_moist.    */
    /**********************************************************************/
    max_infil = (1.0 + b_infilt) * max_moist_tot;
    if (b_infilt == -1.0) {
        tmp = max_infil;
    }
    else {
        ratio = 1.0 - (moist) / (max_moist_tot);
        if (ratio > 1.0) {
            log_err("SOIL RATIO GREATER THAN 1.0: moisture %f  "
                    "max_moisture %f -> ratio = %f",
                    moist, max_moist_tot, ratio);
        }
        else {
            if (ratio < 0.0) {
                log_err("SOIL RATIO LESS THAN 0.0: moisture %f   "
                        "max_moisture %f -> ratio = %e",
                        moist, max_moist_tot, ratio);
            }
            else {
                ratio = pow(ratio, (1.0 / (b_infilt + 1.0)));
            }
        }
        tmp = max_infil * (1.0 - ratio);
    }

    /**********************************************************************/
    /* Evaporate at potential rate, i.e., Eq.(10) in Liang's derivation.  */
    /**********************************************************************/

    if (tmp >= max_infil) {
        esoil = Epot;
    }
    else {
        /********************************************************************/
        /*  Compute As. 'As' is % area saturated, '1-As' is % area          */
        /*  that is unsaturated.                                            */
        /********************************************************************/

        ratio = tmp / max_infil;
        ratio = 1.0 - ratio;

        if (ratio > 1.0) {
            log_err("EVAP RATIO GREATER THAN 1.0");
        }
        else {
            if (ratio < 0.0) {
                log_err("EVAP RATIO LESS THAN 0.0");
            }
            else {
                if (ratio != 0.0) {
                    ratio = pow(ratio, b_infilt);
                }
            }
        }

        as = 1 - ratio;

        /********************************************************************/
        /*  Compute the beta function in the ARNO evaporation model using   */
        /*  the first 30 terms in the power expansion expression.           */
        /********************************************************************/

        ratio = pow(ratio, (1.0 / b_infilt));

        dummy = 1.0;
        for (num_term = 1; num_term <= 30; num_term++) {
            tmpsum = ratio;
            for (i = 1; i < num_term; i++) {
                tmpsum *= ratio;
            }
            dummy += b_infilt * tmpsum / (b_infilt + num_term);
        }

        beta_asp = as + (1.0 - as) * (1.0 - ratio) * dummy;
        esoil = Epot * beta_asp;
    }

    /***********************************************************************/
    /*  Evaporation cannot exceed available soil moisture.                 */
    /*  Evaporation second soil layer = 0.0                                */
    /***********************************************************************/

    /* only consider positive evaporation; we won't put limits on condensation */
    if (esoil > 0.0) {
        if (moist > resid_moist_tot) {
            /* there is liquid moisture available; cap esoil at available liquid moisture */
            if (esoil > moist - resid_moist_tot) {
                esoil = moist - resid_moist_tot;
            }
        }
        else {
            /* no moisture available; cap esoil at 0 */
            esoil = 0.0;
        }
    }

    tmpesoil = esoil;
    for (lidx = 0; lidx < Nlayers_esoil; lidx++) {
        if (layer[lidx].moist - resid_moist[lidx] > tmpesoil) {
            layer[lidx].esoil = tmpesoil;
            tmpesoil = 0;
        }
        else {
            layer[lidx].esoil = layer[lidx].moist - resid_moist[lidx];
            tmpesoil -= layer[lidx].moist - resid_moist[lidx];
        }
    }
    Evap += esoil / MM_PER_M / delta_t;

    return(Evap);
}
