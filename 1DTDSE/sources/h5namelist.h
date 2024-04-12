/**
 * @file h5namelist.h
 * @brief Contains names used to specify h5 paths
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef H5NAMELIST_H
#define H5NAMELIST_H

extern const char * CTDSE_inputs;
extern const char * CUPRAD_outputs;
extern const char * global_outputs;

extern char * CUPRAD_outputs_Efield;
extern char * CUPRAD_outputs_tgrid;
extern char * CUPRAD_outputs_rgrid;
extern char * CUPRAD_outputs_zgrid;

extern char * CUPRAD_outputs_kz_step;
extern char * CUPRAD_outputs_kr_step;
extern char * CUPRAD_outputs_Nz_max;
extern char * CUPRAD_outputs_Nr_max;


void Init_h5_paths();

#endif