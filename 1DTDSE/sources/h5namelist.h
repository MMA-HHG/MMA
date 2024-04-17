/**
 * @file h5namelist.h
 * @brief Define constants' names used to specify h5 paths
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef H5NAMELIST_H
#define H5NAMELIST_H

// paths to the respective groups
#define CTDSE_INPUTS        "CTDSE/inputs/"
#define CUPRAD_OUTPUTS      "CUPRAD/outputs/"
#define GLOBAL_INPUTS       "global_inputs/"


// automatic definitions
#define CUPRAD_OUTPUTS_EFIELD   CUPRAD_OUTPUTS "output_field"
#define CUPRAD_OUTPUTS_TGRID    CUPRAD_OUTPUTS "tgrid"
#define CUPRAD_OUTPUTS_RGRID    CUPRAD_OUTPUTS "rgrid"
#define CUPRAD_OUTPUTS_ZGRID    CUPRAD_OUTPUTS "zgrid"

#define CTDSE_INPUTS_KZ_STEP    CTDSE_INPUTS "kz_step"
#define CTDSE_INPUTS_KR_STEP    CTDSE_INPUTS "kr_step"
#define CTDSE_INPUTS_NZ_MAX     CTDSE_INPUTS "Nz_max"
#define CTDSE_INPUTS_NR_MAX     CTDSE_INPUTS "Nr_step"

#endif



// extern const char * CTDSE_inputs;
// extern const char * CUPRAD_outputs;
// extern const char * global_outputs;

// extern char * CUPRAD_outputs_Efield;
// extern char * CUPRAD_outputs_tgrid;
// extern char * CUPRAD_outputs_rgrid;
// extern char * CUPRAD_outputs_zgrid;

// extern char * CTDSE_inputs_kz_step;
// extern char * CTDSE_inputs_kr_step;
// extern char * CTDSE_inputs_Nz_max;
// extern char * CTDSE_inputs_Nr_max;


// void Init_h5_paths();
