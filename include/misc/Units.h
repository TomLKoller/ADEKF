//This file defines macros to use different  units with better visibility in the code
// Created by tomlucas on 04.06.20.
//

#ifndef CYLINDEREXAMPLE_SI_UNITS_H
#define CYLINDEREXAMPLE_SI_UNITS_H


//time units base is seconds
#define unit_s 1.
#define unit_hr 3600.
#define unit_min 60.

#define _s
#define _hr *unit_hr
#define _p_s
#define _p_hr /unit_hr
#define _rhr *sqrt(unit_hr)
#define _p_rhr /sqrt(unit_hr)

#define _min *unit_min
#define _p_min /unit_min


//angle units base is rad
#define unit_rad 1.
#define unit_deg M_PI/180.

#define _rad
#define _deg *unit_deg


//distance unit base is m
#define unit_m 1.
#define unit_km 1000.

#define _m
#define _p_m

#define _km *unit_km
#define _p_km /unit_km

#endif //CYLINDEREXAMPLE_SI_UNITS_H
