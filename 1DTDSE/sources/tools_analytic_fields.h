/**
 * @file tools_analytic_fields.h
 * @brief Header containing functions for analytical and numerical fields.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TOOLS_ANALYTIC_FIELDS_H
#define TOOLS_ANALYTIC_FIELDS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "structures.h"

double Afieldflattop1(double, double, double , double , double , double , double , double );
double Afieldflattop1ch(double, double, double , double , double , double , double , double , double , double);
double smootherstep(double , double, double);
double clamp(double, double, double);
double Afieldsin2(double, double, double, double, double, double, double);
double dAfieldsin2(double, double, double, double, double, double, double);
double Primsin2cos(double , double , double , double , double );
double AfieldEsin2(double , double , double , double , double , double , double );
//double EField(double,double,double,double,int,double,double);
double AField(Efield_var, double);
double dAField(Efield_var, double);

#endif