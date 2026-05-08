#include "config.h"
#include "mysound.h"
#include <stdio.h>

int main(void)
{
    CapConfig cfg = {.device = "sysdefault:CARD=Set"};
    CaptureCtx ctx = {0};

    // Initialize mixer handle and element
    open_capture_elem(&ctx, &cfg);
    if (!ctx.mixerHandle || !ctx.mixerElem)
    {
        printf("Failed to open ALSA mixer or find 'Mic' element.\n");
        return 1;
    }

    long value0 = -1;
    if (get_mic_amplification(&ctx, &value0) == 0)
    {
        printf("Current mic amplification: %ld\n", value0);
    }
    {
        const long TWENTY = 20;
        long value = TWENTY;
        if (get_mic_amplification(&ctx, &value) != 0)
        {
            printf("Failed to get mic amplification.\n");
            return 2;
        }

        // Loop to increase in steps of 5 until max
        {
            const long STEP_SIZE = 5;
            long new_value = value + STEP_SIZE;
            if (set_mic_amplification(&ctx, new_value) == 0)
            {
                printf("Set mic amplification to: %ld\n", new_value);
            }
            else
            {
                printf("Failed to set mic amplification to: %ld\n", new_value);
            }
            long check_value = -1;
            if (get_mic_amplification(&ctx, &check_value) == 0)
            {
                printf("Mic amplification after set: %ld\n", check_value);
            }
            else
            {
                printf("Failed to get mic amplification after set.\n");
            }
            value = check_value;
        }

        // Clean up mixer handle
        if (ctx.mixerHandle)
        {
            snd_mixer_close(ctx.mixerHandle);
            ctx.mixerHandle = NULL;
            ctx.mixerElem = NULL;
        }

        return 0;
    }
}
