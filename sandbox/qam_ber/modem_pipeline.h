#ifndef MODEM_PIPELINE_H
#define MODEM_PIPELINE_H

#include <complex.h>
#include "liquid.h"

typedef struct {
    unsigned int   dbps;
    unsigned int   M_data;
    modemcf        mod;

    unsigned int   k;
    unsigned int   m;
    float          beta;
    firinterp_crcf interp;
    firdecim_crcf  decim;
} modem_pipeline_t;

int  modem_pipeline_create(modem_pipeline_t *mp, modulation_scheme ms,
                           unsigned int k, unsigned int m, float beta);
void modem_pipeline_destroy(modem_pipeline_t *mp);
void modem_pipeline_reset(modem_pipeline_t *mp);

/* TX: modulate + RRC interpolate. Writes k samples to tx_buf. */
void modem_pipeline_tx(modem_pipeline_t *mp, unsigned int data,
                       float complex *tx_buf);

/* RX: RRC decimate + signal conditioning. Returns one symbol-rate sample. */
void modem_pipeline_rx(modem_pipeline_t *mp, float complex *rx_buf,
                       float complex *y_out);

/* Decode: demodulate received symbol-rate samples. */
void modem_pipeline_decode(modem_pipeline_t *mp,
                           const float complex *rx_samp,
                           unsigned int num_symbols, unsigned int *dec_out);

/* Training: feed one RX symbol + known reference to adapt RX conditioning.
   ref = known TX constellation point for this symbol.
   Call after pipeline settles (i.e. after total_delay symbols). */
void modem_pipeline_train_step(modem_pipeline_t *mp,
                               float complex *rx_buf,
                               float complex ref);

/* Finalize training: freeze RX conditioning for measurement phase. */
void modem_pipeline_train_done(modem_pipeline_t *mp);

/* Generate N deterministic training TX symbols into ref[N]. */
void modem_pipeline_gen_train_ref(modem_pipeline_t *mp,
                                  float complex *ref, unsigned int N);

/* Total pipeline delay in symbols. */
unsigned int modem_pipeline_delay(const modem_pipeline_t *mp);

#endif
