# tico

Audio capture and terminal visualization tool that captures from a microphone and displays a time-domain graph that scrolls across your terminal.

## Features

- Real-time audio capture from microphone
- FFT-based frequency analysis
- Customizable logging options (e.g., average sound per tick with peak detection)
- Terminal-based visualization
- Multi-compiler support (GCC, Clang 15, Clang 19)
- Debug and Release build configurations

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
| `debug-clang15` | Clang 15 | Debug | Debug build with Clang 15 |
| `release-clang15` | Clang 15 | Release | Release build with Clang 15 |
| `debug-clang19` | Clang 19 | Debug | Debug build with Clang 19 |
| `release-clang19` | Clang 19 | Release | Release build with Clang 19 |

#### Build Examples

```bash
# GCC builds
cmake --preset debug-gcc
cmake --build --preset debug-gcc

cmake --preset release-gcc
cmake --build --preset release-gcc

# Clang 15 builds
cmake --preset debug-clang15
cmake --build --preset debug-clang15

cmake --preset release-clang15
cmake --build --preset release-clang15

# Clang 19 builds
cmake --preset debug-clang19
cmake --build --preset debug-clang19

cmake --preset release-clang19
cmake --build --preset release-clang19
```

### Workflow Presets

For convenience, workflow presets combine configure → build → test steps:

```bash
# Run full workflow (configure + build + test)
cmake --workflow --preset debug
cmake --workflow --preset release
cmake --workflow --preset debug-gcc
cmake --workflow --preset release-gcc
cmake --workflow --preset debug-clang15
cmake --workflow --preset release-clang15
cmake --workflow --preset debug-clang19
cmake --workflow --preset release-clang19
```

## Testing

Run tests using CTest presets:

```bash
# Test with GCC
ctest --preset debug_ctest
ctest --preset release_ctest

# Test with Clang 15
ctest --preset debug_ctest-clang15
ctest --preset release_ctest-clang15

# Test with Clang 19
ctest --preset debug_ctest-clang19
ctest --preset release_ctest-clang19

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
  - `build/clang15-debug/` - Clang 15 debug
  - `build/clang15-release/` - Clang 15 release
  - `build/clang19-debug/` - Clang 19 debug
  - `build/clang19-release/` - Clang 19 release

### Compiler Versions

- GCC (default system version)
- Clang 15
- Clang 19

The project requires C17 standard support and exports compile commands for IDE integration.

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

(Add your license information here)