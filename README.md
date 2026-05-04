# tico

Audio capture and terminal visualization tool that captures from a microphone and displays a time-domain graph that scrolls across your terminal.

## Features

- Real-time audio capture from microphone (robust automatic device selection)
- FFT-based frequency analysis
- Customizable logging options (e.g., average sound per tick with peak detection)
- Terminal-based visualization
- Multi-compiler support (GCC, latest available Clang)
- Debug and Release build configurations

## Microphone Device Selection

By default, tico automatically selects the most appropriate microphone device using robust heuristics:

- **USB microphones are preferred**: If a sysdefault device with a USB description is found, it is used automatically.
- **Sysdefault fallback**: If no USB device is found, any sysdefault device is used.
- **Final fallback**: If no sysdefault device is found, the standard "default" ALSA device is used.

This logic ensures that USB microphones are prioritized for best compatibility and quality, but the program will work out-of-the-box on most Linux systems with ALSA.

You can override the device selection by specifying the device name with the `-d` option.

## Prerequisites

Before building, install the required dependencies:

```bash
apt-get install libfftw3-dev libasound2-dev cmake
```

- **libfftw3-dev**: FFTW3 library for FFT calculations
- **libasound2-dev**: ALSA sound library for audio capture
- **cmake**: Build system (3.25.0 or later)

## Building

This project uses CMake with preset configurations for different compilers and build types.

### Quick Start (GCC, Debug)

```bash
cmake --preset debug
cmake --build --preset debug
```

### Available Presets

#### Configure Presets

| Preset | Compiler | Build Type | Description |
|--------|----------|-----------|-------------|
| `debug` | GCC | Debug | Debug build with default GCC |
| `release` | GCC | Release | Release build with default GCC |
| `debug-gcc` | GCC | Debug | Explicit GCC debug build |
| `release-gcc` | GCC | Release | Explicit GCC release build |
| `debug-clang` | Clang (auto-detected latest) | Debug | Debug build with newest installed Clang |
| `release-clang` | Clang (auto-detected latest) | Release | Release build with newest installed Clang |
| `doxygen` | Clang (auto-detected latest) | Debug | Generate API documentation with Doxygen |

#### Build Examples

```bash
# GCC builds
cmake --preset debug-gcc
cmake --build --preset debug-gcc

cmake --preset release-gcc
cmake --build --preset release-gcc

# Clang latest builds
cmake --preset debug-clang
cmake --build --preset debug-clang

cmake --preset release-clang
cmake --build --preset release-clang
```

### Workflow Presets

For convenience, workflow presets combine configure → build → test steps:

```bash
# Run full workflow (configure + build + test)
cmake --workflow --preset debug
cmake --workflow --preset release
cmake --workflow --preset debug-gcc
cmake --workflow --preset release-gcc
cmake --workflow --preset debug-clang
cmake --workflow --preset release-clang
cmake --workflow --preset doxygen
```

### Quality Workflow

Run formatting checks, clang-tidy, clang build-and-test, and a GCC build in one command:

```bash
cmake --workflow --preset quality
```

The `quality` workflow runs these steps:

- Configure (`debug-clang`)
- Format check (`format-check`)
- Static analysis (`clang-tidy`)
- Build (`debug-clang`)
- Test (`debug_ctest-clang`)
- GCC build (`quality-gcc-build` -> runs `debug-gcc` configure + build)
- Documentation (`doxygen`)

## Documentation

API documentation is generated using [Doxygen](https://www.doxygen.nl/).

### Doxygen Documentation Workflows

You can generate API documentation as part of a workflow or as a standalone step:

- **Debug/Clang docs:**
  ```bash
  cmake --workflow --preset doxygen
  ```
  Output: `build/clang-debug/html/`

- **Release/GCC docs:**
  ```bash
  cmake --workflow --preset release
  ```
  (Runs doxygen as the last step; output: `build/release/html/`)

- **Release/Clang docs:**
  ```bash
  cmake --workflow --preset release-clang
  ```
  (Runs doxygen as the last step; output: `build/clang-release/html/`)

You can also use the build presets directly:

- `cmake --build --preset doxygen` (debug/clang)
- `cmake --build --preset doxygen-release` (release/gcc)
- `cmake --build --preset doxygen-release-clang` (release/clang)

> **Note:** The generated HTML documentation will be found in the `html/` directory inside the build directory for the preset you use.

## Testing

Run tests using CTest presets:

```bash
# Test with GCC
ctest --preset debug_ctest
ctest --preset release_ctest

# Test with Clang latest
ctest --preset debug_ctest-clang
ctest --preset release_ctest-clang

# Test with explicit GCC
ctest --preset debug_ctest-gcc
ctest --preset release_ctest-gcc
```

## Usage

### Basic Capture

```bash
./tico
```

### With Logging Options

```bash
# Average sound of a tick per tooth with peak detection
./tico -j 15 -p peak
```

## Development

### Project Structure

- Build outputs are organized by compiler and build type:
  - `build/debug/` - GCC debug build
  - `build/release/` - GCC release build
  - `build/gcc-debug/` - Explicit GCC debug
  - `build/gcc-release/` - Explicit GCC release
  - `build/clang-debug/` - Clang debug
  - `build/clang-release/` - Clang release

### Compiler Versions

- GCC (default system version)
- Clang (newest available version detected at configure time)

The project requires C11 standard support and exports compile commands for IDE integration.

### Building from Scratch

```bash
# Clean previous builds (optional)
rm -rf build/

# Configure with your chosen preset
cmake --preset debug-gcc

# Build
cmake --build --preset debug-gcc

# Test
ctest --preset debug_ctest-gcc
```

## License

This project is licensed under the GNU General Public License v3.0 (GPLv3).

**Note:** The GPLv3 license is required because this project depends on the FFTW library, which is licensed under GPL. As a result, any software linked with FFTW must also be distributed under the GPL.

See the LICENSE file for the full license text.