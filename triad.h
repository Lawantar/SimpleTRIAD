#ifndef TRIAD_H
#define TRIAD_H

#include <stdint.h>

#define JULDAY_CONST 1721422.5f
#define SEC_PER_YEAR 31556952.0f
#define DEGREES_TO_RADIANS 0.01745329252f

#define SUN_OP_MODE 0x00
#define NADIR_OP_MODE 0x01

#define DIV_ZERO_ERR 0x01
#define NO_SUN_VEC 0x02
#define NO_CORRECT_DATA 0x04
#define UNKNOWN_OP_MODE 0x08
#define NO_MAG_DATA 0x10

#define NUM_OF_OM 6 //просто рандмное значение для компиля :)

/* IGNORE */
typedef struct {
    uint32_t* cmd; // +00, [ComReg] 4 байта. Команда ориентации.
    uint8_t* orient_status; // +04, [StatReg] 4 байта. Статус работы ориеентации. Глобальный: остановлена, вработе, ошибка, завершена.
    uint8_t* control; // +08, [ComReg] 128 байт. Данные для команды (зависит от режима
    uint8_t* mag_func_status; // +12, [StatReg] 128 байт. Статус работы ориеентации. Локальный, для функции.
    uint16_t data_err_flag; // +16 Флаги наличия ошибок в данных СОП. Бит выставлен - в данных ошибка.
    uint16_t res0; // +18, резерв
    uint8_t alg_stop_flag; // +20, Флаг завершения алгоритма ориентации.
    uint8_t alg_err_flag; // +21, Флаг ошибки алгоритма ориентации.
    uint16_t res1; // +22, резерв
    uint32_t call_time; // +24, Момент времени текущего вызова функции
    uint32_t last_call_time; // +28, Момент времени предыдущего вызова функции
    uint32_t last_exit_time; // +32, Момент времени предыдущего выхода из функции
    uint32_t delay_ms; // +36, Запрос задержки перед следующим вызовом
    void (*Func)(void); // +40, сама функция, вызываемая в данный момент
} mag_orientation_ctrl_t;


typedef struct {
    uint32_t alt; // Координаты
    uint32_t lat;
    uint32_t lon;
    uint32_t GNSS_ts;
    uint16_t GNSS_status;
} sens_pos_t;

typedef struct {
    float zenith[NUM_OF_OM]; // Зенит и
    float azimuth[NUM_OF_OM]; // Азимут Солнца
    float centers[NUM_OF_OM][2]; // Положение центра пятна по x и y.
    float Sun[3]; // Вектор на Солнце (расчитаный старым алгоритмом, усреднением углов)
    uint32_t timestamp;
    uint16_t status;
    uint8_t req;
    uint8_t ready;
} sens_OM_SS_t;

typedef struct {
    float vectors[NUM_OF_OM][3]; // Вектор на Землю
    float Nadir[3]; //TODO Нету в СОПе
    uint32_t timestamp;
    uint16_t status;
    uint8_t req;
    uint8_t ready;
} sens_OM_HS_t;

typedef struct {
    float mag_Gs[NUM_OF_OM][3]; // Магнитометры ДСГ
    float gyro_dps[NUM_OF_OM][3]; // ДУСы ДСГ
    float accel_G[NUM_OF_OM][3]; // Акселерометры ДСГ, G (~9.81 mpss)
    uint32_t timestamp;
    uint8_t status_mag[NUM_OF_OM]; // Статус: == 0 => ошибок нет
    uint8_t status_GA[NUM_OF_OM]; // Статус: == 0 => ошибок нет
    uint8_t req; // Запрос данных
    uint8_t ready; // Выставляется при чтении данных из сенсора
} sens_OM_GAM_t;

typedef struct {
    float data_Gs[3]; // Магнитометр СОП
    uint32_t timestamp;
    uint16_t status; // Статус: == 0 => ошибок нет
    uint8_t req; // Запрос данных
    uint8_t ready; // Выставляется при чтении данных из сенсора
} sens_ADCS_mag_t;
typedef struct {
    float gyro_dps[3]; // ДУС СОП
    float accel_G[3]; // Акселерометр СОП, G (~9.81 mpss)
    uint32_t timestamp;
    uint16_t status;
    uint8_t req;
    uint8_t ready;
} sens_ADCS_GA_t;


typedef struct {
    sens_ADCS_mag_t sens_ADCS_mag; // +00
    sens_ADCS_GA_t sens_ADCS_GA; // +20
    sens_OM_GAM_t sens_OM_GAM; // +52
    sens_OM_SS_t sens_OM_SS; // +286
    sens_OM_HS_t sens_OM_HS; // +402
    sens_pos_t sens_pos; // +482
} sensorics_data_t; // [764 bytes]

typedef struct {
    float mag_moment[3];
    uint32_t coil_time;
} mag_control_t;
/* IGNORE */

typedef struct {
    float gcrf2itrf_matrix[9];
    float triad_matrix[9];
    float b_model[3];
    float sun_model[3];
    float nadir_model[3];
    float jul_date;
    uint8_t op_mode;
} model_t;

void ConvertDateToJD(void);

void GetSunDirection(void);

void GetNadirFromGeo(void);

void GCRFToITRF(void);

void normalize(float v[3]);

void cross(const float a[3], const float b[3], float result[3]);

float dot(const float a[3], const float b[3]);

void vecMatInv(const float matrix[9], const float vector[3]);

void matInv(const float matrix[9], float inv[9]);

uint8_t simpleTriad(void);

void prisma(void);

#endif