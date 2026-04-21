#ifndef CHANNEL_SIM_H
#define CHANNEL_SIM_H

#include <complex.h>
#include "liquid.h"

typedef struct {
    firfilt_crcf bw_filt;
    unsigned int bw_filt_len;
    firfilt_crcf mp_filt;
    unsigned int mp_len;
    unsigned int k;
} channel_sim_t;

int  channel_sim_create(channel_sim_t *ch, float sym_rate, unsigned int k,
                        float bw_cutoff);
void channel_sim_destroy(channel_sim_t *ch);
void channel_sim_reset(channel_sim_t *ch);

void channel_sim_apply(channel_sim_t *ch, const float complex *in,
                       unsigned int n, float complex *out);
void channel_sim_add_noise(float complex *buf, unsigned int n, float nstd);

unsigned int channel_sim_delay(const channel_sim_t *ch);

#endif
