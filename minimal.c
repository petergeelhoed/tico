#include <alsa/asoundlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
int read_16000_samples_i16_nonblock(snd_pcm_t* cap, int16_t* out)
{
    const unsigned TARGET = 16000;
    unsigned collected = 0;

    while (collected < TARGET)
    {
        snd_pcm_sframes_t got =
            snd_pcm_readi(cap, out + collected, TARGET - collected);

        if (got == -EAGAIN)
        {
            /* No data available right now */
            /* Yield or sleep very briefly */
            usleep(1000); /* 1 ms */
            continue;
        }

        if (got == -EPIPE || got == -ESTRPIPE)
        {
            /* XRUN or suspend */
            if (snd_pcm_recover(cap, (int)got, 1) < 0)
                return -1;
            continue;
        }

        if (got < 0)
        {
            fprintf(stderr, "ALSA read failed: %s\n", snd_strerror((int)got));
            return -1;
        }

        collected += (unsigned)got;
    }

    return (int)TARGET;
}

int main(void)
{
    snd_pcm_t* cap = NULL;
    snd_pcm_hw_params_t* hw = NULL;
    unsigned rate = 48000;
    int dir = 0;

    /* Open device in NON-BLOCKING mode */
    if (snd_pcm_open(&cap,
                     "default",
                     SND_PCM_STREAM_CAPTURE,
                     SND_PCM_NONBLOCK) < 0)
    {
        fprintf(stderr, "Cannot open capture device\n");
        return 1;
    }

    snd_pcm_hw_params_malloc(&hw);
    snd_pcm_hw_params_any(cap, hw);

    snd_pcm_hw_params_set_access(cap, hw, SND_PCM_ACCESS_RW_INTERLEAVED);
    snd_pcm_hw_params_set_format(cap, hw, SND_PCM_FORMAT_S16_LE);
    snd_pcm_hw_params_set_channels(cap, hw, 1);
    snd_pcm_hw_params_set_rate_near(cap, hw, &rate, &dir);

    if (snd_pcm_hw_params(cap, hw) < 0)
    {
        fprintf(stderr, "Cannot set HW params\n");
        return 1;
    }

    snd_pcm_hw_params_free(hw);
    snd_pcm_prepare(cap);

    int16_t samples[16000];

    printf("Recording (non-blocking)...\n");

    if (read_16000_samples_i16_nonblock(cap, samples) != 16000)
    {
        fprintf(stderr, "Capture failed\n");
        return 1;
    }

    printf("Done.\n");

    for (int i = 0; i < 16000; ++i)
        printf("%d\n", samples[i]);

    snd_pcm_close(cap);
    return 0;
}
