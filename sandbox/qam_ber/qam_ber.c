#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <getopt.h>
#include "liquid.h"
#include "modem_pipeline.h"
#include "channel_sim.h"

static void run_training(modem_pipeline_t *mp, channel_sim_t *ch,
                         unsigned int N_train, float SNRdB)
{
    unsigned int k = mp->k;
    unsigned int total_delay = modem_pipeline_delay(mp) + channel_sim_delay(ch);
    float nstd = powf(10.0f, -SNRdB / 20.0f) * sqrtf((float)k);

    float complex *ref = malloc(N_train * sizeof(float complex));
    modem_pipeline_gen_train_ref(mp, ref, N_train);

    modem_pipeline_reset(mp);
    channel_sim_reset(ch);
    modem_pipeline_train_begin(mp, N_train);

    float complex tx_buf[k], ch_buf[k];
    for (unsigned int i = 0; i < N_train + total_delay; i++) {
        unsigned int d = (i < N_train) ? ((i * 7 + 3) % mp->M_data) : 0;
        modem_pipeline_tx(mp, d, tx_buf);
        channel_sim_apply(ch, tx_buf, k, ch_buf);
        channel_sim_add_noise(ch_buf, k, nstd);

        if (i >= total_delay && (i - total_delay) < N_train)
            modem_pipeline_train_step(mp, ch_buf, ref[i - total_delay]);
        else {
            float complex dummy;
            modem_pipeline_rx(mp, ch_buf, &dummy);
        }
    }
    modem_pipeline_train_done(mp);
    free(ref);
}

int main(int argc, char *argv[])
{
    unsigned int k           = 8, m = 7;
    float        beta        = 0.3f;
    float        sym_rate    = 3000.0f;
    float        bw_cutoff   = 1900.0f;
    unsigned int num_symbols = 5000000;
    unsigned int N_train     = 8000;
    float        SNRdB_min   = 28.75f;
    float        SNRdB_max   = 28.75f;
    float        SNRdB_step  = 0.25f;
    const char  *mod_str     = "qam256";

    int dopt;
    while ((dopt = getopt(argc, argv, "hm:n:s:x:d:c:t:")) != EOF) {
        switch (dopt) {
        case 'h':
            printf("qam_ber - BER with BW-limited + multipath channel\n");
            printf("  -m <scheme> modulation scheme [qam256]\n");
            printf("  -n <num>    symbols [5000000]\n");
            printf("  -t <num>    training symbols [8000]\n");
            printf("  -s/-x/-d    SNR min/max/step dB [28.75/28.75/0.25]\n");
            printf("  -c <Hz>     channel BW cutoff [1900]\n");
            return 0;
        case 'm': mod_str     = optarg; break;
        case 'n': num_symbols = atoi(optarg); break;
        case 't': N_train     = atoi(optarg); break;
        case 's': SNRdB_min   = atof(optarg); break;
        case 'x': SNRdB_max   = atof(optarg); break;
        case 'd': SNRdB_step  = atof(optarg); break;
        case 'c': bw_cutoff   = atof(optarg); break;
        default:  return 1;
        }
    }

    modulation_scheme ms = liquid_getopt_str2mod(mod_str);
    if (ms == LIQUID_MODEM_UNKNOWN) {
        fprintf(stderr, "error: unknown modulation scheme '%s'\n", mod_str);
        return 1;
    }

    modem_pipeline_t mp;
    modem_pipeline_create(&mp, ms, k, m, beta);

    channel_sim_t ch;
    channel_sim_create(&ch, sym_rate, k, bw_cutoff);

    unsigned int total_delay = modem_pipeline_delay(&mp) + channel_sim_delay(&ch);

    unsigned int  *tx_data = malloc(num_symbols * sizeof(unsigned int));
    float complex *rx_samp = malloc(num_symbols * sizeof(float complex));
    unsigned int  *dec_out = malloc(num_symbols * sizeof(unsigned int));

    printf("# %s\n", mod_str);
    printf("# channel           : BW %.0f Hz + multipath [1, +.03, -.025, +.02, -.015]\n", bw_cutoff);
    float sim_duration = (float)num_symbols / sym_rate;
    printf("# sim duration      : %.3f s (%u symbols @ %.0f sym/s)\n",
           sim_duration, num_symbols, sym_rate);
    printf("# effective bitrate : %.0f bps (%.1f kbps)\n",
           sym_rate * mp.dbps, sym_rate * mp.dbps / 1000.0);
    printf("# %8s %8s %10s %12s\n", "SNR [dB]", "errors", "bits", "BER");

    float complex tx_buf[k], ch_buf[k];

    for (float SNRdB = SNRdB_min; SNRdB <= SNRdB_max + 0.001f; SNRdB += SNRdB_step) {
        run_training(&mp, &ch, N_train, SNRdB);

        for (unsigned int i = 0; i < num_symbols; i++)
            tx_data[i] = rand() % mp.M_data;

        modem_pipeline_reset(&mp);
        channel_sim_reset(&ch);
        float nstd = powf(10.0f, -SNRdB / 20.0f) * sqrtf((float)k);
        unsigned int rx_n = 0;

        for (unsigned int i = 0; i < num_symbols + total_delay; i++) {
            unsigned int d = (i < num_symbols) ? tx_data[i] : 0;
            modem_pipeline_tx(&mp, d, tx_buf);
            channel_sim_apply(&ch, tx_buf, k, ch_buf);
            channel_sim_add_noise(ch_buf, k, nstd);
            float complex y;
            modem_pipeline_rx(&mp, ch_buf, &y);
            if (i >= total_delay && rx_n < num_symbols)
                rx_samp[rx_n++] = y;
        }

        modem_pipeline_decode(&mp, rx_samp, num_symbols, dec_out);

        unsigned long errs = 0, bits = 0;
        for (unsigned int i = 0; i < num_symbols; i++) {
            errs += count_bit_errors(tx_data[i], dec_out[i]);
            bits += mp.dbps;
        }
        printf("  %8.2f %8lu %10lu  %12.4e\n",
               SNRdB, errs, bits, (float)errs / (float)bits);
        if (errs == 0) break;
    }

    free(tx_data); free(rx_samp); free(dec_out);
    modem_pipeline_destroy(&mp);
    channel_sim_destroy(&ch);
    return 0;
}
