#ifndef JULIALIB_PAINTMIX_H
#define JULIALIB_PAINTMIX_H
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct JLWStatus {
    int32_t code;
    uint8_t message[256];
} JLWStatus;
typedef struct ModelInfo {
    JLWStatus status;
    int32_t abi_version;
    int32_t format_version;
    int32_t grid_n;
    int32_t channel_count;
    uint32_t flags;
    uint32_t forward_crc32;
    uint32_t inverse_crc32;
    uint64_t model_id_lo;
    uint64_t model_id_hi;
} ModelInfo;
typedef struct CVector_borrowed_Float64 {
    int64_t dims[1];
    double* data;
} CVector_borrowed_Float64;
typedef struct CVector_borrowed_UInt32 {
    int64_t dims[1];
    uint32_t* data;
} CVector_borrowed_UInt32;
typedef struct CVector_borrowed_Float32 {
    int64_t dims[1];
    float* data;
} CVector_borrowed_Float32;

JLWStatus paintmix_encode_f32(CVector_borrowed_Float32 rgb, CVector_borrowed_Float32 latent);
JLWStatus paintmix_weighted_mix_f64(CVector_borrowed_Float64 colors, CVector_borrowed_Float64 weights, CVector_borrowed_Float64 out, int64_t count);
JLWStatus paintmix_mix_f32(CVector_borrowed_Float32 a, CVector_borrowed_Float32 b, float t, CVector_borrowed_Float32 out);
JLWStatus paintmix_decode_f32(CVector_borrowed_Float32 latent, CVector_borrowed_Float32 rgb);
JLWStatus paintmix_bulk_mix_f64(CVector_borrowed_Float64 a, CVector_borrowed_Float64 b, CVector_borrowed_Float64 t, CVector_borrowed_Float64 out, int64_t count);
ModelInfo paintmix_model_info();
int32_t paintmix_abi_version();
JLWStatus paintmix_encode_f64(CVector_borrowed_Float64 rgb, CVector_borrowed_Float64 latent);
JLWStatus paintmix_table_crc32(int32_t index, CVector_borrowed_UInt32 out);
JLWStatus paintmix_weighted_mix_f32(CVector_borrowed_Float32 colors, CVector_borrowed_Float32 weights, CVector_borrowed_Float32 out, int64_t count);
JLWStatus paintmix_bulk_mix_f32(CVector_borrowed_Float32 a, CVector_borrowed_Float32 b, CVector_borrowed_Float32 t, CVector_borrowed_Float32 out, int64_t count);
JLWStatus paintmix_mix_f64(CVector_borrowed_Float64 a, CVector_borrowed_Float64 b, double t, CVector_borrowed_Float64 out);
JLWStatus paintmix_decode_f64(CVector_borrowed_Float64 latent, CVector_borrowed_Float64 rgb);
#endif // JULIALIB_PAINTMIX_H
