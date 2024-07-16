#include "triad.h"
#include <math.h>

// Возвращает юлианскую дату.
// Аргументы:
// secs - Кол-во секунд, прошедших от 01.01.2000

F_64 ConvertDateToJD(UI_32 secs) {
    F_64 year = 2000.0 + (F_64) secs / SEC_PER_YEAR;
    F_64 cor = 2 - (year / 100) + (year / 400);
    return JULDAY_CONST + ((year - 1) * 365.25) + cor;
}

// Рассчитывает направление на Солнце, т.е. единичный вектор от Земли к Солнцу. Результат приведен в GCRF.
// Аргументы:
// julDate - Юлианская дата
// sun_dir - Массив из 3х элементов для сохранения направления на Солнце

void GetSunDirection(F_64 julDate, F_64 *sun_dir) {
    F_64 rad = M_PI / 180.0;
    F_64 t_cent = (julDate - 2451545.0) / 36525.0;
    F_64 lambda_m = 280.4606184 + 36000.77005361 * t_cent;
    F_64 M = (357.5277233 + 35999.05034 * t_cent) * rad;
    F_64 lambda_elliptic = (lambda_m + 1.914666471 * sin(M) + 0.019994643 * sin(2 * M)) * rad;
    F_64 eps = (23.439291 - 0.0130042 * t_cent) * rad;
    sun_dir[0] = cos(lambda_elliptic);
    sun_dir[1] = cos(eps) * sin(lambda_elliptic);
    sun_dir[2] = sin(eps) * sin(lambda_elliptic);
}

// Рассчитывает магнитное поле в ITRF
// Аргументы:
// r - Радиус-вектор спутника в ITRF [м]
// earth_dipole - Магнитуда земного диполя [м^3*кг/с^2/А]
// mag_field_GF - Массив из 3х элементов для сохранения магнитного поля

void getMagFieldInclinedDipole(const F_64 *r, F_64 earth_dipole, F_64 *mag_field_GF) {
    F_64 lambda0 = 107.32 * M_PI / 180.0;
    F_64 delta0 = 9.41 * M_PI / 180.0;
    F_64 k[3];
    k[0] = cos(lambda0) * sin(delta0);
    k[1] = sin(lambda0) * sin(delta0);
    k[2] = -cos(delta0);
    F_64 norm_r = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    F_64 k_dot_r = k[0] * r[0] + k[1] * r[1] + k[2] * r[2];
    for (SI_32 i = 0; i < 3; i++) {
        mag_field_GF[i] = -earth_dipole * (k[i] * norm_r * norm_r - 3.0 * k_dot_r * r[i]) / pow(norm_r, 5);
    }
}

// Вычисляет вектор направления в центр Зелми из текущего положения спутника по данным ГЛОНАСС
// Аргументы:
// lat - Широта [град]
// lon - Долгота [град]
// height - Высота по данным ГЛОНАСС [м]
// *nadir - Массив из 3х элементов для сохранения вектора в надир [м]

void GetNadirFromGeo(F_64 lat, F_64 lon, F_64 height, F_64 *nadir) {
    F_64 clat = cos(lat * DEGREES_TO_RADIANS);
    F_64 slat = sin(lat * DEGREES_TO_RADIANS);
    F_64 clon = cos(lon * DEGREES_TO_RADIANS);
    F_64 slon = sin(lon * DEGREES_TO_RADIANS);
    F_64 N = 6378137.0 / sqrt(1.0 - 6.6943799901377997e-3 * slat * slat);
    nadir[0] = -1 * (N + height) * clat * clon;
    nadir[1] = -1 * (N + height) * clat * slon;
    nadir[2] = -1 * (N * (1.0 - 6.6943799901377997e-3) + height) * slat;
}

// Вычисляет матрицу поворота из инерциальной системы координат (GCRF) в систему координат Гринвича (ITRF)
// Аргументы:
// matrix - Массив из 9 элементов для сохранения матрицы поворота
// julDate - Юлианская дата

void GCRFToITRF(F_64 *matrix, F_64 julDate) {
    const F_64 DJ00 = 2451545.0;
    F_64 t_jd_from_J2000 = julDate - DJ00;
    F_64 theta = 2 * M_PI * (0.7790572732640 + 1.00273781191135448 * t_jd_from_J2000);
    matrix[0] = cos(theta);
    matrix[1] = sin(theta);
    matrix[2] = 0.0;
    matrix[3] = -sin(theta);
    matrix[4] = cos(theta);
    matrix[5] = 0.0;
    matrix[6] = 0.0;
    matrix[7] = 0.0;
    matrix[8] = 1.0;
}

void normalize(F_64 v[3]) {
    F_64 n = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n != 0) {
        for (SI_32 i = 0; i < 3; i++) {
            v[i] /= n;
        }
    }
}

void cross(const F_64 a[3], const F_64 b[3], F_64 result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

F_64 dot(const F_64 a[3], const F_64 b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void matInv(const F_64 matrix[9], F_64 inv[9]) {
    F_64 det = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
               matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
               matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
    if (det != 0) {
        F_64 invDet = 1.0 / det;
        inv[0] = invDet * (matrix[4] * matrix[8] - matrix[5] * matrix[7]);
        inv[1] = invDet * (matrix[2] * matrix[7] - matrix[1] * matrix[8]);
        inv[2] = invDet * (matrix[1] * matrix[5] - matrix[2] * matrix[4]);
        inv[3] = invDet * (matrix[5] * matrix[6] - matrix[3] * matrix[8]);
        inv[4] = invDet * (matrix[0] * matrix[8] - matrix[2] * matrix[6]);
        inv[5] = invDet * (matrix[3] * matrix[2] - matrix[0] * matrix[5]);
        inv[6] = invDet * (matrix[3] * matrix[7] - matrix[6] * matrix[4]);
        inv[7] = invDet * (matrix[6] * matrix[1] - matrix[0] * matrix[7]);
        inv[8] = invDet * (matrix[0] * matrix[4] - matrix[3] * matrix[1]);
    } else {
        //Error handle
    }
}

// Реализует алгоритм TRIAD, который рассчитывает матрицу ориентации КА
// Аргументы:
// S_meas - Вектор измерения магнитометра (массив 3 элемента)
// B_meas - Вектор измерения солнечного датчика (массив 3 элемента)
// S_model - Вектор направления на солнце, рассчитаный по моделе (массив 3 элемента)
// B_model - Вектор геомагнитного поля, рассчитаный по моделе (массив 3 элемента)
// A - Матрица ориентации КА (массив 9 элементов)

void simpleTriad(const F_64 S_meas[3], F_64 B_meas[3], const F_64 S_model[3], F_64 B_model[3], F_64 *matrix) {
    normalize(B_meas);
    normalize(B_model);
    F_64 cross_S_meas_B_meas[3];
    F_64 cross_S_model_B_model[3];
    cross(S_meas, B_meas, cross_S_meas_B_meas);
    cross(S_model, B_model, cross_S_model_B_model);
    F_64 D_meas[9] = {S_meas[0], B_meas[0], cross_S_meas_B_meas[0],
                      S_meas[1], B_meas[1], cross_S_meas_B_meas[1],
                      S_meas[2], B_meas[2], cross_S_meas_B_meas[2]};
    F_64 D_model[9] = {S_model[0], B_model[0], cross_S_model_B_model[0],
                       S_model[1], B_model[1], cross_S_model_B_model[1],
                       S_model[2], B_model[2], cross_S_model_B_model[2]};
    F_64 D_model_inv[9];
    matInv(D_model, D_model_inv);
    for (SI_32 i = 0; i < 3; i++) {
        for (SI_32 j = 0; j < 3; j++) {
            matrix[(i * 3) + j] = 0.0;
            for (SI_32 k = 0; k < 3; k++) {
                matrix[(i * 3) + j] += D_meas[(i * 3) + k] * D_model_inv[(k * 3) + j];
            }
        }
    }
    F_64 e1[3] = {matrix[0], matrix[3], matrix[6]};
    F_64 e2[3] = {matrix[1], matrix[4], matrix[7]};
    F_64 e3[3] = {matrix[2], matrix[5], matrix[8]};
    F_64 dot_e1_e2 = dot(e1, e2);
    for (SI_32 i = 0; i < 3; i++) {
        e2[i] -= e1[i] * dot_e1_e2;
    }
    normalize(e2);
    F_64 dot_e1_e3 = dot(e1, e3);
    F_64 dot_e2_e3 = dot(e2, e3);
    for (SI_32 i = 0; i < 3; i++) {
        e3[i] -= e1[i] * dot_e1_e3 + e2[i] * dot_e2_e3;
    }
    normalize(e3);
    for (SI_32 i = 0; i < 3; i++) {
        matrix[i * 3] = e1[i];
        matrix[(i * 3) + 1] = e2[i];
        matrix[(i * 3) + 2] = e3[i];
    }
}

