#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "modem_pipeline.h"

int modem_pipeline_create(modem_pipeline_t *mp, modulation_scheme ms,
                          unsigned int k, unsigned int m, float beta)
{
    memset(mp, 0, sizeof(*mp));
    mp->mod = modemcf_create(ms);
    mp->dbps = modemcf_get_bps(mp->mod);
    mp->M_data = 1u << mp->dbps;
    mp->k = k;  mp->m = m;  mp->beta = beta;

    unsigned int h_len = 2 * k * m + 1;
    float h[h_len];
    liquid_firdes_rrcos(k, m, beta, 0, h);
    mp->interp = firinterp_crcf_create(k, h, h_len);

    float g[h_len];
    for (unsigned int i = 0; i < h_len; i++) g[i] = h[h_len - 1 - i];
    mp->decim = firdecim_crcf_create(k, g, h_len);
    firdecim_crcf_set_scale(mp->decim, 1.0f / (float)k);
    return 0;
}

void modem_pipeline_destroy(modem_pipeline_t *mp)
{
    if (mp->mod)    modemcf_destroy(mp->mod);
    if (mp->interp) firinterp_crcf_destroy(mp->interp);
    if (mp->decim)  firdecim_crcf_destroy(mp->decim);
}

void modem_pipeline_reset(modem_pipeline_t *mp)
{
    firinterp_crcf_reset(mp->interp);
    firdecim_crcf_reset(mp->decim);
}

unsigned int modem_pipeline_delay(const modem_pipeline_t *mp)
{
    return 2 * mp->m;
}

void modem_pipeline_tx(modem_pipeline_t *mp, unsigned int data, float complex *tx_buf)
{
    float complex x;
    modemcf_modulate(mp->mod, data, &x);
    firinterp_crcf_execute(mp->interp, x, tx_buf);
}

void modem_pipeline_rx(modem_pipeline_t *mp, float complex *rx_buf,
                       float complex *y_out)
{
    firdecim_crcf_execute(mp->decim, rx_buf, y_out);
}

void modem_pipeline_train_begin(modem_pipeline_t *mp, unsigned int N_total)
{
    (void)mp;
    (void)N_total;
}

void modem_pipeline_train_step(modem_pipeline_t *mp,
                               float complex *rx_buf,
                               float complex ref)
{
    (void)ref;
    float complex dummy;
    firdecim_crcf_execute(mp->decim, rx_buf, &dummy);
}

void modem_pipeline_train_done(modem_pipeline_t *mp)
{
    (void)mp;
}

void modem_pipeline_decode(modem_pipeline_t *mp, const float complex *rx_samp,
                           unsigned int num_symbols, unsigned int *dec_out)
{
    modemcf dem = modemcf_copy(mp->mod);
    for (unsigned int t = 0; t < num_symbols; t++)
        modemcf_demodulate(dem, rx_samp[t], &dec_out[t]);
    modemcf_destroy(dem);
}

void modem_pipeline_gen_train_ref(modem_pipeline_t *mp, float complex *ref, unsigned int N)
{
    for (unsigned int i = 0; i < N; i++) {
        unsigned int d = (i * 7 + 3) % mp->M_data;
        modemcf_modulate(mp->mod, d, &ref[i]);
    }
}
