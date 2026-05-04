#pragma once

// Returns the first 'plughw:' device from arecord -L output (skipping comments/blank lines)
const char* get_default_device(void);
#include "config.h"
#include "myarr.h"
#include <alsa/asoundlib.h> // IWYU pragma: export
#include <stdio.h>

/** * @brief Context for audio capture, containing ALSA parameters and
 * configuration.
 */
typedef struct CaptureCtx
{
    /* ALSA */
    snd_pcm_t* cap;
    unsigned int rate;
    snd_pcm_uframes_t periodSize;
    snd_pcm_uframes_t bufferSize;
    size_t ArrayLength; // frames per processing block
} CaptureCtx;
/** * @brief Initializes the ALSA audio capture device with the specified
 * format, device name, and sample rate.
 *
 * @param format The audio format to use (e.g., SND_PCM_FORMAT_S16_LE).
 * @param device The name of the ALSA device to open (e.g., "default").
 * @param rate A pointer to an unsigned int where the actual sample rate will
 * be stored. This may be modified if the requested rate is not supported.
 * @return A pointer to the initialized snd_pcm_t capture handle, or NULL on
 * failure.
 */
// NOLINTNEXTLINE(misc-include-cleaner)
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate);

/** * @brief Initializes the audio source for capture based on the provided
 * configuration. This function sets up the ALSA capture device and retrieves
 * the actual sample rate being used.
 *
 * @param cfg A pointer to the CapConfig structure containing the desired
 * configuration for audio capture.
 * @param actualRate A pointer to an unsigned int where the actual sample rate
 * will be stored after initialization.
 * @return An integer status code (0 for success, non-zero for failure).
 */
int initAudioSource(CapConfig* cfg, unsigned int* actualRate);

/** * @brief Reads audio data from the capture device or a file into a buffer.
 *
 * This function reads audio samples from the ALSA capture device or from a
 * specified input file, depending on the context. The read samples are stored
 * in the provided buffer.
 *
 * @param derivative A pointer to an integer array where the read audio samples
 * will be stored.
 * @param ArrayLength The number of audio samples to read into the buffer.
 * @param fpInput A pointer to a FILE object representing the input file to read
 * from. If NULL, audio will be read from the capture device.
 * @param ctx A pointer to the CaptureCtx structure containing the capture
 * context and configuration.
 * @param buffer16 A pointer to an int16_t array where the raw audio samples
 * will be stored before being processed into the derivative array.
 * @return An integer status code (0 for success, non-zero for failure).
 */
int readBufferOrFile(int* derivative,
                     size_t ArrayLength,
                     FILE* fpInput,
                     CaptureCtx* ctx,
                     int16_t* buffer16);

/** * @brief Retrieves audio data from the capture context or input file and
 * processes it into a derivative array.
 *
 * This function reads audio samples from the ALSA capture device or from a
 * specified input file, depending on the context. The read samples are stored
 * in the provided output buffer, and a derivative array is computed based on
 * the raw audio samples.
 *
 * @param fpInput A pointer to a FILE object representing the input file to read
 * from. If NULL, audio will be read from the capture device.
 * @param derivative A pointer to a myarr structure where the computed
 * derivative values will be stored.
 * @param ctx A pointer to the CaptureCtx structure containing the capture
 * context and configuration.
 * @param out A pointer to an int16_t array where the raw audio samples will be
 * stored before being processed into the derivative array.
 * @return An integer status code (0 for success, non-zero for failure).
 */
int getData(FILE* fpInput,
            struct myarr* derivative,
            CaptureCtx* ctx,
            int16_t* out);

/** * @brief Cleans up and releases resources associated with the audio capture
 * context.
 *
 * This function resets the CaptureCtx structure, effectively releasing any
 * resources associated with the audio capture context. It is important to call
 * this function when the capture context is no longer needed to prevent
 * resource leaks.
 *
 * @param ctx A pointer to the CaptureCtx structure to be cleaned up.
 */
void captureTeardown(CaptureCtx* ctx);

/** * @brief Sets up the ALSA audio capture device based on the provided context
 * and configuration.
 *
 * This function configures the ALSA capture device according to the parameters
 * specified in the CaptureCtx structure and the CapConfig configuration. It
 * initializes the capture device and prepares it for audio data capture.
 *
 * @param ctx A pointer to the CaptureCtx structure containing the capture
 * context and configuration.
 * @param cfg A pointer to the CapConfig structure containing the desired
 * configuration for audio capture.
 * @param rate The sample rate to be used for audio capture.
 * @return An integer status code (0 for success, non-zero for failure).
 */
int captureSetup(CaptureCtx* ctx, CapConfig* cfg, unsigned int rate);

/** * @brief Reads audio samples from the ALSA capture device into a buffer.
 *
 * This function reads a specified number of audio samples from the ALSA
 * capture device and stores them in the provided output buffer. It handles
 * potential errors such as EAGAIN (no data available) and XRUN (buffer
 * underrun) by implementing appropriate recovery strategies.
 *
 * @param cap A pointer to the snd_pcm_t capture handle representing the ALSA
 * capture device.
 * @param ArrayLength The number of audio samples to read into the output
 * buffer.
 * @param out A pointer to an int16_t array where the read audio samples will be
 * stored.
 * @return The number of samples successfully read, or a negative error code on
 * failure.
 */
int readSamples(snd_pcm_t* cap, size_t ArrayLength, int16_t* out);
