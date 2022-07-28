#include <iostream>
#include <memory>
#include <vector>
#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <math.h>
#include <time.h>

static constexpr int WINDOW_WIDTH = 640;
static constexpr int WINDOW_HEIGHT = 480;

static const uint32_t   g_sample_rate = 44100;
static float g_sample_time = 1 / (float)g_sample_rate;
float g_time = 0;
static SDL_AudioSpec    g_sdl_audio_spec;
static SDL_AudioSpec    g_sdl_obtained_audio_spec;

enum Waveform {
    SINE,
    TRIANGLE,
    SQUARE,
    SAWTOOTH,
    NOISE
};

Waveform g_wf;
float g_freq;

class Tone
{
public:
    enum ToneName
    {
        C = 0,
        CS,
        D,
        DS,
        E,
        F,
        FS,
        G,
        GS,
        A,
        AS,
        H,
        B = ToneName::H,
    };

    /// Initializes table of frequencies
    Tone();

    /// @param[in] octave Octave from 0..9
    /// @param[in] tone   Counting from 0: C, C#, D, D#, E. F, F#, G, G#, A, A#, H
    float getFrequency(uint8_t octave, ToneName tone);

private:
    std::vector <float> frequencies;
};

Tone::Tone()
{
    double twelfthroot_of_two;
    twelfthroot_of_two = pow(2.0, 1.0 / 12.0);
    frequencies.clear();
    frequencies.resize(10*12);

    double f = 440.0;
    // calculate from C-0 to G#-4
    for (uint8_t i=0; i<(48+9); i++) {
        f /= twelfthroot_of_two;
        frequencies.at(56 - i) = (float)f;
    }

    frequencies.at(57) = 440.0;

    // calculate from A#-4 to H-9
    f = 440.0;
    for (uint8_t i=0; i<(60+2); i++) {
        f *= twelfthroot_of_two;
        frequencies.at(58 + i) = (float)f;
    }

    std::cerr << "C-0: " << getFrequency(0, C) << std::endl;
    std::cerr << "C-3: " << getFrequency(3, C) << std::endl;
    std::cerr << "C-4: " << getFrequency(4, C) << std::endl;
    std::cerr << "A-4: " << getFrequency(4, A) << std::endl;
    std::cerr << "C-5: " << getFrequency(5, C) << std::endl;
    std::cerr << "A-5: " << getFrequency(5, A) << std::endl;
    std::cerr << "B-5: " << getFrequency(5, B) << std::endl;
    std::cerr << "C-6: " << getFrequency(6, C) << std::endl;
    std::cerr << "A-6: " << getFrequency(6, A) << std::endl;
    std::cerr << "C-9: " << getFrequency(9, C) << std::endl;
    std::cerr << "H-9: " << getFrequency(9, H) << std::endl;
}

float Tone::getFrequency(uint8_t octave, ToneName tone)
{
    size_t index = octave * 12 + tone;
    if (index >= frequencies.size()) return 1.0f;
    return frequencies.at(index);
}

class Oscillator
{
public:
    Oscillator();

    void setWaveform(Waveform waveform);
    void setFrequency(float frequency);
    float getAmplitude(float time);

    void noteOn(float time);
    void noteOff(float time);

private:
    Waveform waveform;
    float frequency;
    float time_on;
    float time_off;
};

Oscillator::Oscillator()
    : waveform(SINE)
    , frequency(0.0f)
    , time_on(0.0f)
    , time_off(0.0f)
{
}

void Oscillator::setWaveform(Waveform waveform)
{
    this->waveform = waveform;
}

void Oscillator::setFrequency(float frequency)
{
    if (frequency > 15.0f) {
        this->frequency = frequency;
    }
}

float Oscillator::getAmplitude(float time)
{
    if (time_on < 0.0001f) return 0.0f;
    if (time_off > 0.0001f && time_off < time) return 0.0f;

    float t = time - time_on;
    float T = 1 / this->frequency;
    float mod_t = fmod(t, T);

    float value = 0.0f;

    switch (this->waveform) {
    case SINE:
        value =  sinf(2 * 3.141f * t * this->frequency);
        break;

    case TRIANGLE:
    {
        if (mod_t <= T/4) value = 4*mod_t/T;
        else if (mod_t <= 3*T/4) value = 2.0f - 4*mod_t/T;
        else value = -4.0f + 4*mod_t/T;
        break;
    }

    case SQUARE:
    {
        if (mod_t < T/2) value = 1.0f;
        else value = -1.0f;
        break;
    }

    case SAWTOOTH:
    {
        value = -1.0f + 2*(mod_t + T/2)/T;
        if (value > 1.0f) value -= 2.0f;
        break;
    }

    case NOISE:
        value = -1.0f + 2.0f * (float)rand() / (float)RAND_MAX;
        break;

    default:
        break;
    }
    return value;
}

void Oscillator::noteOn(float time)
{
    time_on = time;
    time_off = 0.0f;
}

void Oscillator::noteOff(float time)
{
    time_off = time;
}

class ADSR
{
public:
    ADSR(float a, float d, float s, float r);

    float getAmplitude(float time);
    float getOscillatorValue(float time);
    void noteOn(float time);
    void noteOff(float time);
    bool isOn(float time);

private:
    float attack;
    float decay;
    float sustain;
    float release;

    float time_on;
    float time_off;
    bool is_playing;
};

ADSR::ADSR(float a, float d, float s, float r)
    : attack(a)
    , decay(d)
    , sustain(s)
    , release(r)
    , time_on(0.0f)
    , time_off(0.0f)
{
}

void ADSR::noteOn(float time)
{
    time_on = time;
    time_off = 0.0f;
    is_playing = true;
}

void ADSR::noteOff(float time)
{
    time_off = time_on + attack + decay;
    if (time > time_off) time_off = time;
    is_playing = false;
}

bool ADSR::isOn(float time)
{
    return (time_on > 0.0f && (time_off < 0.001f || (time_off > 0.001 && time < time_off)));
}

float ADSR::getAmplitude(float time)
{
    if (time_on < 0.001) {
        return 0.0f;
    }

    if (time_off < 0.001 || time_off > time) {
        float elapsed = time - time_on;

        if (elapsed <= attack) {
            return elapsed / attack;
        }
        else if (elapsed <= (attack + decay)) {
            elapsed -= attack;
            return sustain + (1.0f - sustain) * (1.0f - elapsed / decay);
        }
        else  return sustain;
    }
    else {
        float elapsed = time - time_off;
        if (elapsed <= release) {
            return sustain * (1.0f - elapsed / release);
        }
        else {
            return -0.001f;
        }
    }
}

Tone g_tone;
std::shared_ptr<Oscillator> g_oscillator = nullptr;
std::shared_ptr<ADSR> g_adsr = nullptr;

float oscillator(float t, float frequency, Waveform waveform)
{
    float value = 0.0f;

    switch (waveform) {
    case SINE:
        value =  sinf(2 * 3.141f * t * frequency);
        break;

    case TRIANGLE:
    {
        float T = 1 / frequency;
        float mod_t = fmod(t, T);

        if (mod_t <= T/4) value = 4*mod_t/T;
        else if (mod_t <= 3*T/4) value = 2.0f - 4*mod_t/T;
        else value = -4.0f + 4*mod_t/T;
        break;
    }

    case SQUARE:
    {
        float T = 1 / frequency;
        float mod_t = fmod(t, T);
        if (mod_t < T/2) value = 1.0f;
        else value = -1.0f;
        break;
    }

    case SAWTOOTH:
    {
        float T = 1 / frequency;
        float mod_t = fmod(t, T);

        value = -1.0f + 2*(mod_t + T/2)/T;
        if (value > 1.0f) value -= 2.0f;
        break;
    }

    case NOISE:
        value = -1.0f + 2.0f * (float)rand() / (float)RAND_MAX;
        break;

    default:
        break;
    }
    return value;
}

void make_noise(float t, float & sample_l, float & sample_r)
{
    //sample_l = sinf(2 * 3.141f * t * 440);
    //sample_r = cosf(2 * 3.141f * t * 439);
    //sample_l = sinf(2 * 3.141f * t * (440 + 1 * sin(2 * 3.141f * 10 * t)));
    //sample_r = cosf(2 * 3.141f * t * 439);
    sample_l = 0.0f;
    sample_r = 0.0f;
    if (g_adsr != nullptr && g_oscillator != nullptr) {
        float v = g_adsr->getAmplitude(t);
        if (v >= 0.0f) {
            sample_l = v * g_oscillator->getAmplitude(t);
            sample_r = sample_l;
        }
        else {
            g_oscillator->noteOff(t);
        }
    }
}

void sound_fill_buffer_s16lsb(uint8_t* streambuf, int bufferlength)
{
    int i, bi, num_samples;
    int16_t l_sample = 0, r_sample = 0;

    // The format of the byte stream is signed 16-bit samples in little-endian 
    // byte order. Stereo samples are stored in a LRLRLR ordering.
    // sample     |             1             |             2             |...
    // channel    |    left     |    right    |    left     |    right    |...
    // byte order | low  | high | low  | high | low  | high | low  | high |...
    // stream     | [0]  | [1]  | [2]  | [3]  | [4]  | [5]  | [6]  | [7]  |...

    // Hence num_samples must be divided by four.

    num_samples = bufferlength >> 2;
    bi = 0; // start with byte index at the beginning

    float fl, fr;

    for (i = 0; i<num_samples; i++) {
        make_noise(g_time, fl, fr);
        g_time += g_sample_time;

        l_sample = (int16_t)(fl * 30000);
        r_sample = (int16_t)(fr * 30000);

        streambuf[bi++] = l_sample & 0x00FF;
        streambuf[bi++] = l_sample >> 8;
        streambuf[bi++] = r_sample & 0x00FF;
        streambuf[bi++] = r_sample >> 8;
    }
}

void sound_fill_buffer_f32lsb(uint8_t* streambuf, int bufferlength)
{
    int i, bi, num_samples;
    int16_t l_sample = 0, r_sample = 0;

    // The format of the byte stream is floating point 32-bit samples in little-endian 
    // byte order. Stereo samples are stored in a LRLRLR ordering.
    // sample     |                           1                           |                           2                           |...
    // channel    |            left           |           right           |            left           |           right           |...
    // byte order | low  | mlow | mhig | high | low  | mlow | mhig | high | low  | mlow | mhig | high | low  | mlow | mhig | high |...
    // stream     | [0]  | [1]  | [2]  | [3]  | [4]  | [5]  | [6]  | [7]  | [8]  | [9]  | [A]  | [B]  | [C]  | [D]  | [E]  | [F]  |...

    // Hence num_samples must be divided by eight.

    num_samples = bufferlength >> 3;
    bi = 0; // start with byte index at the beginning

    union
    {
        float f;
        uint32_t u32;
    } left, right;

    for (i = 0; i<num_samples; i++) {
        make_noise(g_time, left.f, right.f);
        g_time += g_sample_time;

        streambuf[bi++] = (left.u32 & 0x000000FF);
        streambuf[bi++] = (left.u32 & 0x0000FF00) >> 8;
        streambuf[bi++] = (left.u32 & 0x00FF0000) >> 16;
        streambuf[bi++] = (left.u32 & 0xFF000000) >> 24;
        streambuf[bi++] = (right.u32 & 0x000000FF);
        streambuf[bi++] = (right.u32 & 0x0000FF00) >> 8;
        streambuf[bi++] = (right.u32 & 0x00FF0000) >> 16;
        streambuf[bi++] = (right.u32 & 0xFF000000) >> 24;
    }
}

void sound_callback(void* arg, uint8_t* streambuf, int bufferlength)
{
    switch (g_sdl_obtained_audio_spec.format) {
    case AUDIO_S16LSB:
        sound_fill_buffer_s16lsb(streambuf, bufferlength);
        break;
    case AUDIO_F32LSB:
        sound_fill_buffer_f32lsb(streambuf, bufferlength);
        break;
    default:
        // unsupported sample format
        memset(streambuf, 0, bufferlength);
        break;
    }
}

SDL_Window* g_window;
SDL_Renderer* g_renderer;

void logAudioDevices()
{
    int i, count = SDL_GetNumAudioDevices(0);

    for (i = 0; i < count; ++i) {
        SDL_Log("Audio device %d: %s", i, SDL_GetAudioDeviceName(i, 0));
    }
}

const char* SDLAudioFormat2String(SDL_AudioFormat f)
{
    switch (f) {
    case AUDIO_S8: return "signed 8 - bit samples";
    case AUDIO_U8: return "unsigned 8 - bit samples";
    case AUDIO_S16LSB: return "signed 16 - bit samples in little - endian byte order";
    case AUDIO_S16MSB: return "signed 16 - bit samples in big - endian byte order";
    case AUDIO_U16LSB: return "unsigned 16 - bit samples in little - endian byte order";
    case AUDIO_U16MSB: return "unsigned 16 - bit samples in big - endian byte order";
    case AUDIO_S32LSB: return "signed 32 - bit integer samples in little - endian byte order";
    case AUDIO_S32MSB: return "unsigned 32 - bit integer samples in big - endian byte order";
    case AUDIO_F32LSB: return "32 - bit floating point samples in little - endian byte order";
    case AUDIO_F32MSB: return "32 - bit floating point samples in big - endian byte order";

    default:
        SDL_Log("Unknown SDL_AudioFormat %04X", f);
        return "unknown audio format";
    }
}

void logAudioSpec(const char* keyword, SDL_AudioSpec* audio_spec)
{
    if (keyword != nullptr) SDL_Log("%s audio spec", keyword);
    if (audio_spec != nullptr) {
        SDL_Log("    channels: %d", audio_spec->channels);
        SDL_Log("    format:   %s", SDLAudioFormat2String(audio_spec->format));
        SDL_Log("    freq:     %d", audio_spec->freq);
        SDL_Log("    samples:  %d", audio_spec->samples);
        SDL_Log("    size:     %d", audio_spec->size);
    }
    else {
        SDL_Log("NULL");
    }
}

void log_output(void*           userdata,
                int             category,
                SDL_LogPriority priority,
                const char*     message)
{
    puts(message);
}

void init()
{
    do {
        SDL_LogSetOutputFunction(log_output, nullptr); 
        SDL_LogSetAllPriority(SDL_LOG_PRIORITY_VERBOSE);

        if (SDL_Init(SDL_INIT_AUDIO | SDL_INIT_VIDEO) < 0) {
            SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s", SDL_GetError());
            break;
        }

        g_window = SDL_CreateWindow("Soundtest",
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        WINDOW_WIDTH, WINDOW_HEIGHT,
                                        SDL_WINDOW_SHOWN);
        if (g_window == nullptr) {
            std::cerr << "SDL_CreateWindow() failed! " << SDL_GetError() << std::endl;
            break;
        }

        g_renderer = SDL_CreateRenderer(g_window, -1, 0);
        if (g_renderer == nullptr) {
            std::cerr << "SDL_CreateRenderer() failed! " << SDL_GetError() << std::endl;
            break;
        }

        SDL_SetRenderDrawBlendMode(g_renderer, SDL_BLENDMODE_BLEND);

        logAudioDevices();

        g_sdl_audio_spec.freq = g_sample_rate;
        g_sdl_audio_spec.format = AUDIO_S16;
        g_sdl_audio_spec.channels = 2;
        g_sdl_audio_spec.samples = 2048;
        g_sdl_audio_spec.callback = (SDL_AudioCallback)sound_callback;
        g_sdl_audio_spec.userdata = NULL;
        logAudioSpec("wanted", &g_sdl_audio_spec);

        // initialize audio
        if (SDL_OpenAudio(&g_sdl_audio_spec, &g_sdl_obtained_audio_spec) < 0) {
            SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open audio: %s", SDL_GetError());
            break;
        }
        logAudioSpec("obtained", &g_sdl_obtained_audio_spec);
        SDL_PauseAudio(1);

    } while (0);
}

void drawLine(int sx, int sy, int ex, int ey)
{
    SDL_RenderDrawLine(g_renderer, sx, sy, ex, ey);
}

void drawWaveform(Waveform wf)
{
    const float f = g_freq;
    const int N = 100;
    float interval = 1.0f / ((float)N * f);

    int sx = 0, sy = WINDOW_HEIGHT / 2; // axis crosspoint (0/0) of the waveform diagram in the window
    float time = 0.0f;

    int last_x = sx;
    int last_y = sy;
    for (int i=0; i<4*N; i++)
    {
        int x = sx + 160 * i / N;
        int y = (int)(sy - 200 * oscillator(time, f, wf));
        drawLine(last_x, last_y, x, y);
        time += interval;
        last_x = x;
        last_y = y;
    }
}


void noteOn(uint8_t octave, Tone::ToneName tone)
{
    if (g_adsr->isOn(g_time)) return;
    g_freq = g_tone.getFrequency(octave, tone);
    g_oscillator->setFrequency(g_freq);
    g_adsr->noteOn(g_time);
    g_oscillator->noteOn(g_time);
}

void noteOff(uint8_t octave, Tone::ToneName tone)
{
    (void)octave;
    (void)tone;
    g_adsr->noteOff(g_time);
}

#ifdef WINDOWS
int WinMain()
#else
int main(void)
#endif
{
    g_oscillator = std::shared_ptr<Oscillator>(new Oscillator());
    g_adsr = std::shared_ptr<ADSR>(new ADSR(0.1f, 0.2f, 0.4f, 0.4f));
    g_wf = SINE;
    g_freq = 1.0f;

    srand((unsigned int)time(NULL));
    init();

    SDL_Event event;
    bool running = true;

    SDL_PauseAudio(0);
    while (running) {
        SDL_SetRenderDrawColor(g_renderer, 0, 0, 0, 255);
        SDL_RenderClear(g_renderer);
        SDL_SetRenderDrawColor(g_renderer, 255, 255, 255, 255);
        drawWaveform(g_wf);
        SDL_RenderPresent(g_renderer);
        
        while (SDL_WaitEventTimeout(&event, 10)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                case SDLK_1:
                    g_wf = SINE;
                    if (g_oscillator != nullptr) g_oscillator->setWaveform(SINE);
                    break;
                case SDLK_2:
                    g_wf = TRIANGLE;
                    if (g_oscillator != nullptr) g_oscillator->setWaveform(TRIANGLE);
                    break;
                case SDLK_3:
                    g_wf = SQUARE;
                    if (g_oscillator != nullptr) g_oscillator->setWaveform(SQUARE);
                    break;
                case SDLK_4:
                    g_wf = SAWTOOTH;
                    if (g_oscillator != nullptr) g_oscillator->setWaveform(SAWTOOTH);
                    break;
                case SDLK_5:
                    g_wf = NOISE;
                    if (g_oscillator != nullptr) g_oscillator->setWaveform(NOISE);
                    break;

                case SDLK_y: noteOn(4, Tone::ToneName::C); break;
                case SDLK_x: noteOn(4, Tone::ToneName::D); break;
                case SDLK_c: noteOn(4, Tone::ToneName::E); break;
                case SDLK_v: noteOn(4, Tone::ToneName::F); break;
                case SDLK_b: noteOn(4, Tone::ToneName::G); break;
                case SDLK_n: noteOn(4, Tone::ToneName::A); break;
                case SDLK_m: noteOn(4, Tone::ToneName::H); break;
                case SDLK_COMMA: noteOn(5, Tone::ToneName::C); break;

                case SDLK_ESCAPE: running = false; break;

                default: break;
                }
                break;

            case SDL_KEYUP:
                switch (event.key.keysym.sym) {
                case SDLK_y: noteOff(4, Tone::ToneName::C); break;
                case SDLK_x: noteOff(4, Tone::ToneName::D); break;
                case SDLK_c: noteOff(4, Tone::ToneName::E); break;
                case SDLK_v: noteOff(4, Tone::ToneName::F); break;
                case SDLK_b: noteOff(4, Tone::ToneName::G); break;
                case SDLK_n: noteOff(4, Tone::ToneName::A); break;
                case SDLK_m: noteOff(4, Tone::ToneName::H); break;
                case SDLK_COMMA: noteOff(5, Tone::ToneName::C); break;
                default: break;
                }
                break;

            case SDL_QUIT:
                running = false;
                break;

            default:
                break;
            }
        }
    }
    SDL_PauseAudio(1);
    SDL_Quit();
    return 0;
}
