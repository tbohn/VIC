/******************************************************************************
 * @section DESCRIPTION
 *
 * Various functions for manipulating cell storages and fluxes
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
 * @brief    Various functions for manipulating cell storages and fluxes
 *****************************************************************************/

/******************************************************************************
 * @brief    Deep copy of a cell_data_struct.
 *           Assumes destination structure has been allocated already.
 *****************************************************************************/
void
copy_cell_data(cell_data_struct *cell1,
               cell_data_struct *cell2)
{

  extern option_struct options;
  unsigned short l;

  /* Copy all scalars */
  *cell2 = *cell1;
  /* Loop over layers */
  for (l = 0; l < options.Nlayer; l++) {
    /* Copy all scalars */
    (*cell2).layer[l] = (*cell1).layer[l];
  }

}

/******************************************************************************
 * @brief    Deep copy of a veg_var_struct
 *           Assumes destination structure has been allocated already.
 *****************************************************************************/
void
copy_veg_var(veg_var_struct *veg_var1,
             veg_var_struct *veg_var2)
{

  extern option_struct options;
  unsigned short l;

  /* Copy all scalars */
  *veg_var2 = *veg_var1;
  if (options.CARBON) {
    /* Loop over Ncanopy */
    for (l = 0; l < options.Ncanopy; l++) {
      (*veg_var2).NscaleFactor[l] = (*veg_var1).NscaleFactor[l];
      (*veg_var2).aPARLayer[l] = (*veg_var1).aPARLayer[l];
      (*veg_var2).CiLayer[l] = (*veg_var1).CiLayer[l];
      (*veg_var2).rsLayer[l] = (*veg_var1).rsLayer[l];
    }
  }

}

/******************************************************************************
 * @brief    Deep copy of a snow_data_struct.
 *           Assumes destination structure has been allocated already.
 *****************************************************************************/
void
copy_snow_data(snow_data_struct *snow1,
               snow_data_struct *snow2)
{

  extern option_struct options;

  /* Copy all scalars */
  *snow2 = *snow1;

}

/******************************************************************************
 * @brief    Deep copy of an energy_bal_struct.
 *           Assumes destination structure has been allocated already.
 *****************************************************************************/
void
copy_energy_bal(energy_bal_struct *energy1,
                energy_bal_struct *energy2)
{

  extern option_struct options;

  /* Copy all scalars */
  *energy2 = *energy1;

}

/******************************************************************************
 * @brief    Compute weighted average of two cell_data_structs
 *****************************************************************************/
void
wtavg_cell_data(cell_data_struct *cell1,
                double            wt1,
                cell_data_struct *cell2,
                double            wt2,
                bool              do_fluxes,
                cell_data_struct *cell_out)
{

  extern option_struct options;
  unsigned short l;
  unsigned short i;

  // States
  for (i = 0; i < 2; i++) {
    (*cell_out).aero_resist[i] = 1 / ( wt1 / (*cell1).aero_resist[i] + wt2 / (*cell2).aero_resist[i] );
  }
  (*cell_out).asat = wt1 * (*cell1).asat + wt2 * (*cell2).asat;
  (*cell_out).CLitter = wt1 * (*cell1).CLitter + wt2 * (*cell2).CLitter;
  (*cell_out).CInter = wt1 * (*cell1).CInter + wt2 * (*cell2).CInter;
  (*cell_out).CSlow = wt1 * (*cell1).CSlow + wt2 * (*cell2).CSlow;
  (*cell_out).rootmoist = wt1 * (*cell1).rootmoist + wt2 * (*cell2).rootmoist;
  (*cell_out).wetness = wt1 * (*cell1).wetness + wt2 * (*cell2).wetness;
  (*cell_out).zwt = wt1 * (*cell1).zwt + wt2 * (*cell2).zwt;
  (*cell_out).zwt_lumped = wt1 * (*cell1).zwt_lumped + wt2 * (*cell2).zwt_lumped;
  /* Loop over layers */
  for (l = 0; l < options.Nlayer; l++) {
    (*cell_out).layer[l].Cs = wt1 * (*cell1).layer[l].Cs + wt2 * (*cell2).layer[l].Cs;
    (*cell_out).layer[l].T = wt1 * (*cell1).layer[l].T + wt2 * (*cell2).layer[l].T;
    for (i = 0; i < MAX_FROST_AREAS; i++) {
      (*cell_out).layer[l].ice[i] = wt1 * (*cell1).layer[l].ice[i] + wt2 * (*cell2).layer[l].ice[i];
    }
    (*cell_out).layer[l].kappa = wt1 * (*cell1).layer[l].kappa + wt2 * (*cell2).layer[l].kappa;
    (*cell_out).layer[l].moist = wt1 * (*cell1).layer[l].moist + wt2 * (*cell2).layer[l].moist;
    (*cell_out).layer[l].phi = wt1 * (*cell1).layer[l].phi + wt2 * (*cell2).layer[l].phi;
    (*cell_out).layer[l].zwt = wt1 * (*cell1).layer[l].zwt + wt2 * (*cell2).layer[l].zwt;
  }

  if (do_fluxes) {
    // Fluxes
    (*cell_out).baseflow = wt1 * (*cell1).baseflow + wt2 * (*cell2).baseflow;
    (*cell_out).irr_applied = wt1 * (*cell1).irr_applied + wt2 * (*cell2).irr_applied;
    (*cell_out).irr_demand = wt1 * (*cell1).irr_demand + wt2 * (*cell2).irr_demand;
    (*cell_out).irr_run_used = wt1 * (*cell1).irr_run_used + wt2 * (*cell2).irr_run_used;
    (*cell_out).irr_with_used = wt1 * (*cell1).irr_with_used + wt2 * (*cell2).irr_with_used;
    (*cell_out).inflow = wt1 * (*cell1).inflow + wt2 * (*cell2).inflow;
    (*cell_out).pot_evap = wt1 * (*cell1).pot_evap + wt2 * (*cell2).pot_evap;
    (*cell_out).runoff = wt1 * (*cell1).runoff + wt2 * (*cell2).runoff;
    (*cell_out).RhLitter = wt1 * (*cell1).RhLitter + wt2 * (*cell2).RhLitter;
    (*cell_out).RhLitter2Atm = wt1 * (*cell1).RhLitter2Atm + wt2 * (*cell2).RhLitter2Atm;
    (*cell_out).RhInter = wt1 * (*cell1).RhInter + wt2 * (*cell2).RhInter;
    (*cell_out).RhSlow = wt1 * (*cell1).RhSlow + wt2 * (*cell2).RhSlow;
    (*cell_out).RhTot = wt1 * (*cell1).RhTot + wt2 * (*cell2).RhTot;
    /* Loop over layers */
    for (l = 0; l < options.Nlayer; l++) {
      (*cell_out).layer[l].esoil = wt1 * (*cell1).layer[l].esoil + wt2 * (*cell2).layer[l].esoil;
      (*cell_out).layer[l].evap = wt1 * (*cell1).layer[l].evap + wt2 * (*cell2).layer[l].evap;
      (*cell_out).layer[l].transp = wt1 * (*cell1).layer[l].transp + wt2 * (*cell2).layer[l].transp;
    }
  }

}

/******************************************************************************
 * @brief    Compute weighted average of two veg_var_structs
 *****************************************************************************/
void
wtavg_veg_var(veg_var_struct *veg_var1,
              double          wt1,
              veg_var_struct *veg_var2,
              double          wt2,
              bool            do_fluxes,
              veg_var_struct *veg_var_out)
{

  extern option_struct options;
  extern parameters_struct param;
  unsigned short l;

  // States
  (*veg_var_out).albedo = wt1 * (*veg_var1).albedo + wt2 * (*veg_var2).albedo;
  (*veg_var_out).fcrop = wt1 * (*veg_var1).fcrop + wt2 * (*veg_var2).fcrop;
  (*veg_var_out).firr = wt1 * (*veg_var1).firr + wt2 * (*veg_var2).firr;
  (*veg_var_out).LAI = wt1 * (*veg_var1).LAI + wt2 * (*veg_var2).LAI;
  (*veg_var_out).fcanopy = wt1 * (*veg_var1).fcanopy + wt2 * (*veg_var2).fcanopy;
  (*veg_var_out).Wdew = wt1 * (*veg_var1).Wdew + wt2 * (*veg_var2).Wdew;
  (*veg_var_out).Wdmax = (*veg_var_out).LAI * param.VEG_LAI_WATER_FACTOR;
  (*veg_var_out).rc = 1 / ( wt1 / (*veg_var1).rc + wt2 / (*veg_var2).rc );
  if (options.CARBON) {
    (*veg_var_out).Ci = wt1 * (*veg_var1).Ci + wt2 * (*veg_var2).Ci;
    (*veg_var_out).NPPfactor = wt1 * (*veg_var1).NPPfactor + wt2 * (*veg_var2).NPPfactor;
    (*veg_var_out).AnnualNPP = wt1 * (*veg_var1).AnnualNPP + wt2 * (*veg_var2).AnnualNPP;
    (*veg_var_out).AnnualNPPPrev = wt1 * (*veg_var1).AnnualNPPPrev + wt2 * (*veg_var2).AnnualNPPPrev;
    /* Loop over canopy layers */
    for (l = 0; l < options.Ncanopy; l++) {
      (*veg_var_out).NscaleFactor[l] = wt1 * (*veg_var1).NscaleFactor[l] + wt2 * (*veg_var2).NscaleFactor[l];
      (*veg_var_out).CiLayer[l] = wt1 * (*veg_var1).CiLayer[l] + wt2 * (*veg_var2).CiLayer[l];
      (*veg_var_out).rsLayer[l] = 1 / ( wt1 / (*veg_var1).rsLayer[l] + wt2 / (*veg_var2).rsLayer[l] );
    }
  }

  if (do_fluxes) {
    // Fluxes
    (*veg_var_out).canopyevap = wt1 * (*veg_var1).canopyevap + wt2 * (*veg_var2).canopyevap;
    (*veg_var_out).throughfall = wt1 * (*veg_var1).throughfall + wt2 * (*veg_var2).throughfall;
    if (options.CARBON) {
      (*veg_var_out).aPAR = wt1 * (*veg_var1).aPAR + wt2 * (*veg_var2).aPAR;
      (*veg_var_out).GPP = wt1 * (*veg_var1).GPP + wt2 * (*veg_var2).GPP;
      (*veg_var_out).Rphoto = wt1 * (*veg_var1).Rphoto + wt2 * (*veg_var2).Rphoto;
      (*veg_var_out).Rdark = wt1 * (*veg_var1).Rdark + wt2 * (*veg_var2).Rdark;
      (*veg_var_out).Rmaint = wt1 * (*veg_var1).Rmaint + wt2 * (*veg_var2).Rmaint;
      (*veg_var_out).Rgrowth = wt1 * (*veg_var1).Rgrowth + wt2 * (*veg_var2).Rgrowth;
      (*veg_var_out).Raut = wt1 * (*veg_var1).Raut + wt2 * (*veg_var2).Raut;
      (*veg_var_out).NPP = wt1 * (*veg_var1).NPP + wt2 * (*veg_var2).NPP;
      (*veg_var_out).Litterfall = wt1 * (*veg_var1).Litterfall + wt2 * (*veg_var2).Litterfall;
      /* Loop over canopy layers */
      for (l = 0; l < options.Ncanopy; l++) {
        (*veg_var_out).aPARLayer[l] = wt1 * (*veg_var1).aPARLayer[l] + wt2 * (*veg_var2).aPARLayer[l];
      }
    }
  }

}

/******************************************************************************
 * @brief    Compute weighted average of two snow_data_structs
 *****************************************************************************/
void
wtavg_snow_data(snow_data_struct *snow1,
                bool              overstory1,
                double            wt1,
                snow_data_struct *snow2,
                bool              overstory2,
                double            wt2,
                bool              do_fluxes,
                snow_data_struct *snow_out)
{

  extern option_struct options;
  double sum;

  // States
  (*snow_out).albedo = wt1 * (*snow1).albedo + wt2 * (*snow2).albedo;
  (*snow_out).canopy_albedo = 0;
  sum = 0;
  if (overstory1) {
    (*snow_out).canopy_albedo += wt1 * (*snow1).canopy_albedo;
    sum += wt1;
  }
  if (overstory2) {
    (*snow_out).canopy_albedo += wt2 * (*snow2).canopy_albedo;
    sum += wt2;
  }
  if (sum > 0) {
    (*snow_out).canopy_albedo /= sum;
  }
  else {
    (*snow_out).canopy_albedo = (*snow1).canopy_albedo;
  }
  (*snow_out).coldcontent = wt1 * (*snow1).coldcontent + wt2 * (*snow2).coldcontent;
  (*snow_out).coverage = wt1 * (*snow1).coverage + wt2 * (*snow2).coverage;
  (*snow_out).depth = wt1 * (*snow1).depth + wt2 * (*snow2).depth;
  if ((*snow_out).depth > 0) {
    (*snow_out).density = ( wt1 * (*snow1).density * (*snow1).depth + wt2 * (*snow2).density * (*snow2).depth ) / (*snow_out).depth;
  }
  else {
    (*snow_out).density = 0;
  }
  (*snow_out).last_snow = ( (*snow1).last_snow > (*snow2).last_snow ? (*snow1).last_snow : (*snow2).last_snow );
  (*snow_out).max_snow_depth = wt1 * (*snow1).max_snow_depth + wt2 * (*snow2).max_snow_depth;
  if ((*snow1).MELTING || (*snow2).MELTING) {
    (*snow_out).MELTING = true;
  }
  (*snow_out).pack_water = wt1 * (*snow1).pack_water + wt2 * (*snow2).pack_water;
  if ((*snow1).snow || (*snow2).snow) {
    (*snow_out).snow = true;
  }
  (*snow_out).snow_canopy = wt1 * (*snow1).snow_canopy + wt2 * (*snow2).snow_canopy;
  (*snow_out).store_coverage = wt1 * (*snow1).store_coverage + wt2 * (*snow2).store_coverage;
  if ((*snow1).store_snow || (*snow2).store_snow) {
    (*snow_out).store_snow = true;
  }
  (*snow_out).store_swq = wt1 * (*snow1).store_swq + wt2 * (*snow2).store_swq;
  if ((*snow1).surf_temp_fbflag || (*snow2).surf_temp_fbflag) {
    (*snow_out).surf_temp_fbcount = ( (*snow1).surf_temp_fbcount > (*snow2).surf_temp_fbcount ? (*snow1).surf_temp_fbcount : (*snow2).surf_temp_fbcount );
    (*snow_out).surf_temp_fbflag = 1;
  }
  (*snow_out).surf_water = wt1 * (*snow1).surf_water + wt2 * (*snow2).surf_water;
  (*snow_out).swq = wt1 * (*snow1).swq + wt2 * (*snow2).swq;
  (*snow_out).pack_temp = 0;
  (*snow_out).surf_temp = 0;
  (*snow_out).snow_distrib_slope = 0;
  (*snow_out).tmp_int_storage = 0;
  if ((*snow1).swq > 0) {
    (*snow_out).pack_temp += wt1 * (*snow1).pack_temp;
    (*snow_out).surf_temp += wt1 * (*snow1).surf_temp;
    (*snow_out).snow_distrib_slope = wt1 * (*snow1).snow_distrib_slope;
    (*snow_out).tmp_int_storage = wt1 * (*snow1).tmp_int_storage;
  }
  if ((*snow2).swq > 0) {
    (*snow_out).pack_temp += wt2 * (*snow2).pack_temp;
    (*snow_out).surf_temp += wt2 * (*snow2).surf_temp;
    (*snow_out).snow_distrib_slope = wt2 * (*snow2).snow_distrib_slope;
    (*snow_out).tmp_int_storage = wt2 * (*snow2).tmp_int_storage;
  }

  if (do_fluxes) {
    // Fluxes
    (*snow_out).blowing_flux = wt1 * (*snow1).blowing_flux + wt2 * (*snow2).blowing_flux;
    (*snow_out).canopy_vapor_flux = wt1 * (*snow1).canopy_vapor_flux + wt2 * (*snow2).canopy_vapor_flux;
    (*snow_out).mass_error = wt1 * (*snow1).mass_error + wt2 * (*snow2).mass_error;
    (*snow_out).melt = wt1 * (*snow1).melt + wt2 * (*snow2).melt;
    (*snow_out).Qnet = wt1 * (*snow1).Qnet + wt2 * (*snow2).Qnet;
    (*snow_out).surface_flux = wt1 * (*snow1).surface_flux + wt2 * (*snow2).surface_flux;
    (*snow_out).transport = wt1 * (*snow1).transport + wt2 * (*snow2).transport;
    (*snow_out).vapor_flux = wt1 * (*snow1).vapor_flux + wt2 * (*snow2).vapor_flux;
  }

}

/******************************************************************************
 * @brief    Compute weighted average of two energy_bal_structs
 *****************************************************************************/
void
wtavg_energy_bal(energy_bal_struct *energy1,
                 bool               LAKE1,
                 bool               FS_ACTIVE1,
                 double             wt1,
                 energy_bal_struct *energy2,
                 bool               LAKE2,
                 bool               FS_ACTIVE2,
                 double             wt2,
                 double            *Zsum_node,
                 bool               do_fluxes,
                 energy_bal_struct *energy_out)
{

  extern option_struct options;
  unsigned short i;
  double sum;

  // States
  (*energy_out).AlbedoLake = 0;
  sum = 0;
  if (LAKE1) {
    (*energy_out).AlbedoLake += wt1 * (*energy1).AlbedoLake;
    sum += wt1;
  }
  if (LAKE2) {
    (*energy_out).AlbedoLake += wt2 * (*energy2).AlbedoLake;
    sum += wt2;
  }
  if (sum > 0) {
    (*energy_out).AlbedoLake /= sum;
  }
  (*energy_out).AlbedoOver = wt1 * (*energy1).AlbedoOver + wt2 * (*energy2).AlbedoOver;
  (*energy_out).AlbedoUnder = wt1 * (*energy1).AlbedoUnder + wt2 * (*energy2).AlbedoUnder;
  (*energy_out).T1_index = (int) ( wt1 * (*energy1).T1_index + wt2 * (*energy2).T1_index );
  (*energy_out).Tcanopy = wt1 * (*energy1).Tcanopy + wt2 * (*energy2).Tcanopy;
  if ((*energy1).Tcanopy_fbflag || (*energy2).Tcanopy_fbflag) {
    (*energy_out).Tcanopy_fbcount = ( (*energy1).Tcanopy_fbcount > (*energy2).Tcanopy_fbcount ? (*energy1).Tcanopy_fbcount : (*energy2).Tcanopy_fbcount );
    (*energy_out).Tcanopy_fbflag = 1;
  }
  (*energy_out).Tfoliage = wt1 * (*energy1).Tfoliage + wt2 * (*energy2).Tfoliage;
  if ((*energy1).Tfoliage_fbflag || (*energy2).Tfoliage_fbflag) {
    (*energy_out).Tfoliage_fbcount = ( (*energy1).Tfoliage_fbcount > (*energy2).Tfoliage_fbcount ? (*energy1).Tfoliage_fbcount : (*energy2).Tfoliage_fbcount );
    (*energy_out).Tfoliage_fbflag = 1;
  }
  (*energy_out).Tsurf = wt1 * (*energy1).Tsurf + wt2 * (*energy2).Tsurf;
  if ((*energy1).Tsurf_fbflag || (*energy2).Tsurf_fbflag) {
    (*energy_out).Tsurf_fbcount = ( (*energy1).Tsurf_fbcount > (*energy2).Tsurf_fbcount ? (*energy1).Tsurf_fbcount : (*energy2).Tsurf_fbcount );
    (*energy_out).Tsurf_fbflag = 1;
  }
  (*energy_out).unfrozen = wt1 * (*energy1).unfrozen + wt2 * (*energy2).unfrozen;
  for (i = 0; i < 2; i++) {
    (*energy_out).Cs[i] = wt1 * (*energy1).Cs[i] + wt2 * (*energy2).Cs[i];
    (*energy_out).kappa[i] = wt1 * (*energy1).kappa[i] + wt2 * (*energy2).kappa[i];
  }
  for (i = 0; i < MAX_NODES; i++) {
    (*energy_out).Cs_node[i] = wt1 * (*energy1).Cs_node[i] + wt2 * (*energy2).Cs_node[i];
    (*energy_out).ice[i] = wt1 * (*energy1).ice[i] + wt2 * (*energy2).ice[i];
    (*energy_out).kappa_node[i] = wt1 * (*energy1).kappa_node[i] + wt2 * (*energy2).kappa_node[i];
    (*energy_out).moist[i] = wt1 * (*energy1).moist[i] + wt2 * (*energy2).moist[i];
    (*energy_out).T[i] = wt1 * (*energy1).T[i] + wt2 * (*energy2).T[i];
    if ((*energy_out).T[i] < 0) {
      (*energy_out).frozen = true;
    }
    if ((*energy1).T_fbflag[i] || (*energy2).T_fbflag[i]) {
      (*energy_out).T_fbcount[i] = ( (*energy1).T_fbcount[i] > (*energy2).T_fbcount[i] ? (*energy1).T_fbcount[i] : (*energy2).T_fbcount[i] );
      (*energy_out).T_fbflag[i] = 1;
    }
  }
  if(!options.QUICK_FLUX && (FS_ACTIVE1 || FS_ACTIVE2)) {
    find_0_degree_fronts(energy_out, Zsum_node, (*energy_out).T, options.Nnode);
  }

  if (do_fluxes) {
    // Fluxes
    (*energy_out).advected_sensible = wt1 * (*energy1).advected_sensible + wt2 * (*energy2).advected_sensible;
    (*energy_out).advection = wt1 * (*energy1).advection + wt2 * (*energy2).advection;
    (*energy_out).AtmosError = wt1 * (*energy1).AtmosError + wt2 * (*energy2).AtmosError;
    (*energy_out).AtmosLatent = wt1 * (*energy1).AtmosLatent + wt2 * (*energy2).AtmosLatent;
    (*energy_out).AtmosLatentSub = wt1 * (*energy1).AtmosLatentSub + wt2 * (*energy2).AtmosLatentSub;
    (*energy_out).AtmosSensible = wt1 * (*energy1).AtmosSensible + wt2 * (*energy2).AtmosSensible;
    (*energy_out).canopy_advection = wt1 * (*energy1).canopy_advection + wt2 * (*energy2).canopy_advection;
    (*energy_out).canopy_latent = wt1 * (*energy1).canopy_latent + wt2 * (*energy2).canopy_latent;
    (*energy_out).canopy_latent_sub = wt1 * (*energy1).canopy_latent_sub + wt2 * (*energy2).canopy_latent_sub;
    (*energy_out).canopy_refreeze = wt1 * (*energy1).canopy_refreeze + wt2 * (*energy2).canopy_refreeze;
    (*energy_out).canopy_sensible = wt1 * (*energy1).canopy_sensible + wt2 * (*energy2).canopy_sensible;
    (*energy_out).deltaCC = wt1 * (*energy1).deltaCC + wt2 * (*energy2).deltaCC;
    (*energy_out).deltaH = wt1 * (*energy1).deltaH + wt2 * (*energy2).deltaH;
    (*energy_out).error = wt1 * (*energy1).error + wt2 * (*energy2).error;
    (*energy_out).fusion = wt1 * (*energy1).fusion + wt2 * (*energy2).fusion;
    (*energy_out).grnd_flux = wt1 * (*energy1).grnd_flux + wt2 * (*energy2).grnd_flux;
    (*energy_out).latent = wt1 * (*energy1).latent + wt2 * (*energy2).latent;
    (*energy_out).latent_sub = wt1 * (*energy1).latent_sub + wt2 * (*energy2).latent_sub;
    (*energy_out).longwave = wt1 * (*energy1).longwave + wt2 * (*energy2).longwave;
    (*energy_out).LongOverIn = wt1 * (*energy1).LongOverIn + wt2 * (*energy2).LongOverIn;
    (*energy_out).LongUnderIn = wt1 * (*energy1).LongUnderIn + wt2 * (*energy2).LongUnderIn;
    (*energy_out).LongUnderOut = wt1 * (*energy1).LongUnderOut + wt2 * (*energy2).LongUnderOut;
    (*energy_out).melt_energy = wt1 * (*energy1).melt_energy + wt2 * (*energy2).melt_energy;
    (*energy_out).NetLongAtmos = wt1 * (*energy1).NetLongAtmos + wt2 * (*energy2).NetLongAtmos;
    (*energy_out).NetLongOver = wt1 * (*energy1).NetLongOver + wt2 * (*energy2).NetLongOver;
    (*energy_out).NetLongUnder = wt1 * (*energy1).NetLongUnder + wt2 * (*energy2).NetLongUnder;
    (*energy_out).NetShortAtmos = wt1 * (*energy1).NetShortAtmos + wt2 * (*energy2).NetShortAtmos;
    (*energy_out).NetShortOver = wt1 * (*energy1).NetShortOver + wt2 * (*energy2).NetShortOver;
    (*energy_out).NetShortUnder = wt1 * (*energy1).NetShortUnder + wt2 * (*energy2).NetShortUnder;
    (*energy_out).out_long_canopy = wt1 * (*energy1).out_long_canopy + wt2 * (*energy2).out_long_canopy;
    (*energy_out).out_long_surface = wt1 * (*energy1).out_long_surface + wt2 * (*energy2).out_long_surface;
    (*energy_out).refreeze_energy = wt1 * (*energy1).refreeze_energy + wt2 * (*energy2).refreeze_energy;
    (*energy_out).sensible = wt1 * (*energy1).sensible + wt2 * (*energy2).sensible;
    (*energy_out).shortwave = wt1 * (*energy1).shortwave + wt2 * (*energy2).shortwave;
    (*energy_out).ShortOverIn = wt1 * (*energy1).ShortOverIn + wt2 * (*energy2).ShortOverIn;
    (*energy_out).ShortUnderIn = wt1 * (*energy1).ShortUnderIn + wt2 * (*energy2).ShortUnderIn;
    (*energy_out).snow_flux = wt1 * (*energy1).snow_flux + wt2 * (*energy2).snow_flux;
  }

}

