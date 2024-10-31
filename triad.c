#include "triad.h"
#include <math.h>

model_t model;

/* IGNORE */
mag_orientation_ctrl_t mag_orientation_ctrl;
sensorics_data_t sensorics_data;
mag_control_t mag_control;
/* IGNORE */

void ConvertDateToJD(void) {
    float year = 2000.0f + (float) secs / SEC_PER_YEAR; //TODO получать на СОП
    float cor = 2 - (year / 100) + (year / 400);
    model.jul_date = JULDAY_CONST + ((year - 1) * 365.25) + cor;
}

void GetSunDirection(void) {
    float rad = M_PI / 180.0;
    float t_cent = (model.jul_date - 2451545.0f) / 36525.0f;
    float lambda_m = 280.4606184f + 36000.77005361f * t_cent;
    float M = (357.5277233f + 35999.05034f * t_cent) * rad;
    float lambda_elliptic = (lambda_m + 1.914666471f * sinf(M) + 0.019994643f * sinf(2 * M)) * rad;
    float eps = (23.439291f - 0.0130042f * t_cent) * rad;
    model.sun_model[0] = cosf(lambda_elliptic);
    model.sun_model[1] = cosf(eps) * sinf(lambda_elliptic);
    model.sun_model[2] = sinf(eps) * sinf(lambda_elliptic);
}

void GetNadirFromGeo(void) {
    float clat = cosf(sensorics_data.sens_pos.lat * DEGREES_TO_RADIANS);
    float slat = sinf(sensorics_data.sens_pos.lat * DEGREES_TO_RADIANS);
    float clon = cosf(sensorics_data.sens_pos.lon * DEGREES_TO_RADIANS);
    float slon = sinf(sensorics_data.sens_pos.lon * DEGREES_TO_RADIANS);
    float N = 6378137.0f / sqrtf(1.0f - 6.6943799901377997e-3f * slat * slat);
    model.nadir_model[0] = -1 * (N + sensorics_data.sens_pos.alt) * clat * clon;
    model.nadir_model[1] = -1 * (N + sensorics_data.sens_pos.alt) * clat * slon;
    model.nadir_model[2] = -1 * (N * (1.0f - 6.6943799901377997e-3f) + sensorics_data.sens_pos.alt) * slat;
}

void GCRFToITRF(void) {
    const float DJ00 = 2451545.0f;
    float t_jd_from_J2000 = model.jul_date - DJ00;
    float theta = 2.0f * (float) M_PI * (0.7790572732640f + 1.00273781191135448f * t_jd_from_J2000);
    model.gcrf2itrf_matrix[0] = cosf(theta);
    model.gcrf2itrf_matrix[1] = sinf(theta);
    model.gcrf2itrf_matrix[2] = 0.0f;
    model.gcrf2itrf_matrix[3] = -sinf(theta);
    model.gcrf2itrf_matrix[4] = cosf(theta);
    model.gcrf2itrf_matrix[5] = 0.0f;
    model.gcrf2itrf_matrix[6] = 0.0f;
    model.gcrf2itrf_matrix[7] = 0.0f;
    model.gcrf2itrf_matrix[8] = 1.0f;
}

void normalize(float v[3]) {
    float n = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n != 0) {
        for (int8_t i = 0; i < 3; i++) {
            v[i] /= n;
        }
    }
}

void cross(const float a[3], const float b[3], float result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

float dot(const float a[3], const float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void vecMatInv(const float matrix[9], const float vector[3]) {
    float res_vector[3];
    for (int8_t i = 0; i < 3; i++) {
        for (int8_t j = 0; j < 3; j++) {
            res_vector[j] += matrix[(j * 3) + i] * vector[i];
        }
    }
    model.sun_model[0] = res_vector[0];
    model.sun_model[1] = res_vector[1];
    model.sun_model[2] = res_vector[2];
}

void matInv(const float matrix[9], float inv[9]) {
    float det = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
                matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
                matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
    if (det != 0) {
        float invDet = 1.0f / det;
        inv[0] = invDet * (matrix[4] * matrix[8] - matrix[5] * matrix[7]);
        inv[1] = invDet * (matrix[2] * matrix[7] - matrix[1] * matrix[8]);
        inv[2] = invDet * (matrix[1] * matrix[5] - matrix[2] * matrix[4]);
        inv[3] = invDet * (matrix[5] * matrix[6] - matrix[3] * matrix[8]);
        inv[4] = invDet * (matrix[0] * matrix[8] - matrix[2] * matrix[6]);
        inv[5] = invDet * (matrix[3] * matrix[2] - matrix[0] * matrix[5]);
        inv[6] = invDet * (matrix[3] * matrix[7] - matrix[6] * matrix[4]);
        inv[7] = invDet * (matrix[6] * matrix[1] - matrix[0] * matrix[7]);
        inv[8] = invDet * (matrix[0] * matrix[4] - matrix[3] * matrix[1]);
        mag_orientation_ctrl.alg_err_flag &= 0xFF - DIV_ZERO_ERR;
    } else {
        mag_orientation_ctrl.alg_err_flag |= DIV_ZERO_ERR;
        return;
    }
}

uint8_t simpleTriad(void) {
    // Заполнение model_t
    ConvertDateToJD();
    if (sensorics_data.sens_OM_SS.status == 0) {
        GetSunDirection();
        GCRFToITRF();
        vecMatInv(model.gcrf2itrf_matrix, model.sun_model);
        model.op_mode = SUN_OP_MODE;
        mag_orientation_ctrl.alg_err_flag &= 0xFF - NO_CORRECT_DATA - NO_SUN_VEC;
    } else if (sensorics_data.sens_OM_HS.status == 0) {
        GetNadirFromGeo();
        //TODO проверить систему координат для надира
        model.op_mode = NADIR_OP_MODE;
        mag_orientation_ctrl.alg_err_flag &= 0xFF - NO_CORRECT_DATA;
        mag_orientation_ctrl.alg_err_flag |= NO_SUN_VEC;
    } else {
        mag_orientation_ctrl.alg_err_flag |= NO_CORRECT_DATA;
        return 1;
    }

    normalize(sensorics_data.sens_ADCS_mag.data_Gs);
    normalize(model.b_model);

    float cross_V_meas_B_meas[3];
    float cross_V_model_B_model[3];
    float D_meas[9];
    float D_model[9];

    switch (model.op_mode) {
        case SUN_OP_MODE:
            cross(sensorics_data.sens_OM_SS.Sun, sensorics_data.sens_ADCS_mag.data_Gs, cross_V_meas_B_meas);
            cross(model.sun_model, model.b_model, cross_V_model_B_model);
            D_meas[0] = sensorics_data.sens_OM_SS.Sun[0];
            D_meas[1] = sensorics_data.sens_ADCS_mag.data_Gs[0];
            D_meas[2] = cross_V_meas_B_meas[0];
            D_meas[3] = sensorics_data.sens_OM_SS.Sun[1];
            D_meas[4] = sensorics_data.sens_ADCS_mag.data_Gs[1];
            D_meas[5] = cross_V_meas_B_meas[1];
            D_meas[6] = sensorics_data.sens_OM_SS.Sun[2];
            D_meas[7] = sensorics_data.sens_ADCS_mag.data_Gs[2];
            D_meas[8] = cross_V_meas_B_meas[2];

            D_model[0] = model.sun_model[0];
            D_model[1] = model.b_model[0];
            D_model[2] = cross_V_model_B_model[0];
            D_model[3] = model.sun_model[1];
            D_model[4] = model.b_model[1];
            D_model[5] = cross_V_model_B_model[1];
            D_model[6] = model.sun_model[2];
            D_model[7] = model.b_model[2];
            D_model[8] = cross_V_model_B_model[2];
            mag_orientation_ctrl.alg_err_flag &= 0xFF - UNKNOWN_OP_MODE;
            break;

        case NADIR_OP_MODE:
            cross(sensorics_data.sens_OM_HS.Nadir, sensorics_data.sens_ADCS_mag.data_Gs, cross_V_meas_B_meas);
            cross(model.nadir_model, model.b_model, cross_V_model_B_model);
            D_meas[0] = sensorics_data.sens_OM_HS.Nadir[0];
            D_meas[1] = sensorics_data.sens_ADCS_mag.data_Gs[0];
            D_meas[2] = cross_V_meas_B_meas[0];
            D_meas[3] = sensorics_data.sens_OM_HS.Nadir[1];
            D_meas[4] = sensorics_data.sens_ADCS_mag.data_Gs[1];
            D_meas[5] = cross_V_meas_B_meas[1];
            D_meas[6] = sensorics_data.sens_OM_HS.Nadir[2];
            D_meas[7] = sensorics_data.sens_ADCS_mag.data_Gs[2];
            D_meas[8] = cross_V_meas_B_meas[2];

            D_model[0] = model.nadir_model[0];
            D_model[1] = model.b_model[0];
            D_model[2] = cross_V_model_B_model[0];
            D_model[3] = model.nadir_model[1];
            D_model[4] = model.b_model[1];
            D_model[5] = cross_V_model_B_model[1];
            D_model[6] = model.nadir_model[2];
            D_model[7] = model.b_model[2];
            D_model[8] = cross_V_model_B_model[2];
            mag_orientation_ctrl.alg_err_flag &= 0xFF - UNKNOWN_OP_MODE;
            break;

        default:
            mag_orientation_ctrl.alg_err_flag |= UNKNOWN_OP_MODE;
            return 1;
    }

    float D_model_inv[9];
    matInv(D_model, D_model_inv);
    for (int8_t i = 0; i < 3; i++) {
        for (int8_t j = 0; j < 3; j++) {
            model.triad_matrix[(i * 3) + j] = 0.0f;
            for (int8_t k = 0; k < 3; k++) {
                model.triad_matrix[(i * 3) + j] += D_meas[(i * 3) + k] * D_model_inv[(k * 3) + j];
            }
        }
    }
    float e1[3] = {model.triad_matrix[0], model.triad_matrix[3], model.triad_matrix[6]};
    float e2[3] = {model.triad_matrix[1], model.triad_matrix[4], model.triad_matrix[7]};
    float e3[3] = {model.triad_matrix[2], model.triad_matrix[5], model.triad_matrix[8]};
    float dot_e1_e2 = dot(e1, e2);
    for (int8_t i = 0; i < 3; i++) {
        e2[i] -= e1[i] * dot_e1_e2;
    }
    normalize(e2);
    float dot_e1_e3 = dot(e1, e3);
    float dot_e2_e3 = dot(e2, e3);
    for (int8_t i = 0; i < 3; i++) {
        e3[i] -= e1[i] * dot_e1_e3 + e2[i] * dot_e2_e3;
    }
    normalize(e3);
    for (int8_t i = 0; i < 3; i++) {
        model.triad_matrix[i] = e1[i];
        model.triad_matrix[i + 3] = e2[i];
        model.triad_matrix[i + 6] = e3[i];
    }
    return 0;
}

void prisma(void) {
    if (sensorics_data.sens_ADCS_mag.status != 0) {
        mag_orientation_ctrl.alg_err_flag |= NO_MAG_DATA;
        return;
    }
    mag_orientation_ctrl.alg_err_flag &= 0xFF - NO_MAG_DATA;
    if (simpleTriad() == 1) {
        return;
    }

    sensorics_data.sens_ADCS_mag.data_Gs[0] = 1e-6f * sensorics_data.sens_ADCS_mag.data_Gs[0];
    sensorics_data.sens_ADCS_mag.data_Gs[1] = 1e-6f * sensorics_data.sens_ADCS_mag.data_Gs[1];
    sensorics_data.sens_ADCS_mag.data_Gs[2] = 1e-6f * sensorics_data.sens_ADCS_mag.data_Gs[2];
    float moment[3] = {0};
    float e_1[3] = {154.206f, 11.364f, 1};
    float X[3] = {1, 0, 0};
    float Z[3] = {0, 0, 1};

    normalize(e_1);

    float w_0 = 0.1f * (float) M_PI / 180.0f;
    float k = 1e6f;

    float delta[3] = {0, 0, 0};
    delta[0] = e_1[0] - X[0];
    delta[1] = e_1[1] - X[1];
    delta[2] = e_1[2] - X[2];
    normalize(delta);
    delta[0] /= 13.5f;
    delta[1] /= 13.5f;
    delta[2] /= 13.5f;

    float r[3] = {0, 0, 0};

    switch (model.op_mode) {
        case SUN_OP_MODE:
            model.sun_model[0] = model.sun_model[0] + delta[0];
            model.sun_model[1] = model.sun_model[1] + delta[1];
            model.sun_model[2] = model.sun_model[2] + delta[2];
            normalize(model.sun_model);
            cross(model.sun_model, Z, r);
        case NADIR_OP_MODE:
            model.nadir_model[0] = model.nadir_model[0] + delta[0];
            model.nadir_model[1] = model.nadir_model[1] + delta[1];
            model.nadir_model[2] = model.nadir_model[2] + delta[2];
            normalize(model.nadir_model);
            cross(model.nadir_model, Z, r);
        default:
            break;
    }

    float R[3] = {0, 0, 0};
    for (int8_t i = 0; i < 3; i++) {
        R[i] = 0;
        for (int8_t j = 0; j < 3; j++) {
            R[i] += model.triad_matrix[(j * 3) + i] * r[j];
        }
    }

    float w_ref[3];
    for (int8_t i = 0; i < 3; i++) {
        w_ref[i] = w_0 * (e_1[i] + R[i]);
    }

    float w_delta[3];
    for (int8_t i = 0; i < 3; i++) {
        w_delta[i] = sensorics_data.sens_ADCS_GA.gyro_dps[i] - w_ref[i];
    }
    cross(w_delta, sensorics_data.sens_ADCS_mag.data_Gs, moment);
    //TODO Проверка на "norm(moment) > 0.25"
    mag_control.mag_moment[0] = moment[0] * k;
    mag_control.mag_moment[1] = moment[1] * k;
    mag_control.mag_moment[2] = moment[2] * k;
}
