#ifndef TRIAD_H
#define TRIAD_H

#define JULDAY_CONST 1721422.5
#define SEC_PER_YEAR 31556952.0
#define DEGREES_TO_RADIANS 0.01745329252

typedef double F_64;
typedef unsigned long int UI_32;
typedef long int SI_32;

F_64 ConvertDateToJD(UI_32 secs);

void GetSunDirection(F_64 julDate, F_64 *sun_dir);

void getMagFieldInclinedDipole(const F_64 *r, F_64 earth_dipole, F_64 *mag_field_GF);

void GetNadirFromGeo(F_64 lat, F_64 lon, F_64 height, F_64 *nadir);

void GCRFToITRF(F_64* matrix, F_64 t_jd);

void normalize(F_64 v[3]);

void cross(const F_64 a[3], const F_64 b[3], F_64 result[3]);

F_64 dot(const F_64 a[3], const F_64 b[3]);

void matInv(const F_64 matrix[9], F_64 inv[9]);

void simpleTriad(const F_64 S_meas[3], F_64 B_meas[3], const F_64 S_model[3], F_64 B_model[3], F_64 *matrix);

void prisma(F_64 B_bfa[3], F_64 A[9], F_64 S[3], F_64 W[3], F_64 moment[3]);

#endif
