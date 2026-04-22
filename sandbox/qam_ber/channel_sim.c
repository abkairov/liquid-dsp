#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "channel_sim.h"

int channel_sim_create(channel_sim_t *ch, float sym_rate, unsigned int k,
                       float bw_cutoff, unsigned int num_echoes)
{
    memset(ch, 0, sizeof(*ch));
    ch->k = k;

    /* BW-limiting Kaiser LPF */
    float Fs = sym_rate * k;
    ch->bw_filt_len = 16 * k + 1;
    float *h = malloc(ch->bw_filt_len * sizeof(float));
    liquid_firdes_kaiser(ch->bw_filt_len, bw_cutoff / Fs, 50.0f, 0.0f, h);
    float dc = 0;
    for (unsigned int i = 0; i < ch->bw_filt_len; i++) dc += h[i];
    for (unsigned int i = 0; i < ch->bw_filt_len; i++) h[i] /= dc;
    ch->bw_filt = firfilt_crcf_create(h, ch->bw_filt_len);
    free(h);

    /* Multipath: main path + echoes within 4-symbol span */
    if (num_echoes == 0) num_echoes = 4;
    ch->mp_len = 4 * k + 1;
    float *mp = calloc(ch->mp_len, sizeof(float));
    mp[0] = 1.000f;
    for (unsigned int i = 1; i <= num_echoes; i++) {
        unsigned int pos = i * 4 * k / num_echoes;
        float amp = 0.035f - 0.020f * (float)(i - 1) / (float)num_echoes;
        mp[pos] = (i % 2 == 0) ? -amp : amp;
    }
    ch->mp_filt = firfilt_crcf_create(mp, ch->mp_len);
    free(mp);

    return 0;
}

void channel_sim_destroy(channel_sim_t *ch)
{
    if (ch->bw_filt) firfilt_crcf_destroy(ch->bw_filt);
    if (ch->mp_filt) firfilt_crcf_destroy(ch->mp_filt);
}

void channel_sim_reset(channel_sim_t *ch)
{
    firfilt_crcf_reset(ch->bw_filt);
    firfilt_crcf_reset(ch->mp_filt);
}

void channel_sim_apply(channel_sim_t *ch, const float complex *in,
                       unsigned int n, float complex *out)
{
    float complex tmp[n];
    firfilt_crcf_execute_block(ch->bw_filt, (float complex *)in, n, tmp);
    firfilt_crcf_execute_block(ch->mp_filt, tmp, n, out);
}

void channel_sim_add_noise(float complex *buf, unsigned int n, float nstd)
{
    for (unsigned int j = 0; j < n; j++)
        buf[j] += nstd * (randnf() + _Complex_I * randnf()) / sqrtf(2.0f);
}

unsigned int channel_sim_delay(const channel_sim_t *ch)
{
    return (ch->bw_filt_len - 1) / (2 * ch->k);
}
