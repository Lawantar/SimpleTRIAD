#ifndef TRIAD_H
#define TRIAD_H

#define JULDAY_CONST 1721422.5
#define SEC_PER_YEAR 31556952.0

typedef double F_64;
typedef unsigned long int UI_32;
typedef long int SI_32;

F_64 ConvertDateToJD(UI_32 secs);

void GetSunDirection(F_64 julDate, F_64 *sun_dir);

void getMagFieldInclinedDipole(const double *r, double earth_dipole, double *mag_field_GF);

void GCRFToITRF(F_64* matrix, F_64 t_jd);

void normalize(F_64 v[3]);

void cross(const F_64 a[3], const F_64 b[3], F_64 result[3]);

F_64 dot(const F_64 a[3], const F_64 b[3]);

void matInv(const F_64 matrix[9], F_64 inv[9]);

void simpleTriad(const F_64 S_meas[3], F_64 B_meas[3], const F_64 S_model[3], F_64 B_model[3], F_64 *matrix);

#endif
