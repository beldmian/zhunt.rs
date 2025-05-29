# Z-HUNT-3

A high-performance Rust implementation of the Z-DNA formation prediction algorithm, originally developed by Ping-jung Chou under the instruction of Pui S. Ho.

## Overview

Z-HUNT-3 is a computational tool that predicts the formation of Z-DNA (left-handed DNA) in naturally occurring sequences using a thermodynamic approach. This implementation is based on the seminal paper "A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences" published in The EMBO Journal, Vol.5, No.10, pp2737-2744, 1986.

## Features

- **High Performance**: Optimized Rust implementation with parallel processing capabilities
- **Memory Efficient**: Uses dynamic programming to reduce complexity from O(2^n) to O(n)
- **Configurable Threading**: Support for both parallel and sequential processing
- **FASTA Support**: Processes DNA sequences in FASTA format
- **Statistical Analysis**: Provides Z-scores and probability calculations for Z-DNA formation

## Installation

### Prerequisites

- Rust 1.70 or later
- Cargo (comes with Rust)

### Building from Source

```bash
git clone <repository-url>
cd zhunt
cargo build --release
```

The optimized binary will be available at `target/release/zhunt`.

## Usage

```bash
zhunt <window_size> <min_size> <max_size> <filename> [OPTIONS]
```

### Arguments

- `window_size`: Window size in dinucleotides
- `min_size`: Minimum size for analysis
- `max_size`: Maximum size for analysis  
- `filename`: Input sequence file (FASTA format)

### Options

- `-t, --threads <THREADS>`: Number of threads to use (0 = auto-detect) [default: 0]
- `-s, --sequential`: Use sequential processing instead of parallel
- `-h, --help`: Print help information
- `-V, --version`: Print version information

### Examples

```bash
# Basic usage with auto-detected threads
./zhunt 50 6 50 sequence.fasta

# Use specific number of threads
./zhunt 50 6 50 sequence.fasta --threads 8

# Force sequential processing
./zhunt 50 6 50 sequence.fasta --sequential
```

## Output

The program generates a `.Z-SCORE` file containing:

- Position information (start and end positions)
- Sequence length
- Delta linking number
- Slope values
- Probability scores
- DNA sequence segments
- Anti-syn configuration strings

### Output Format

```
filename sequence_length min_size max_size
position end_position length delta_linking slope probability sequence antisyn_config
```

## Algorithm Details

### Z-DNA Formation Prediction

The algorithm uses thermodynamic calculations to predict Z-DNA formation by:

1. **Energy Calculation**: Computing delta BZ energy for dinucleotide pairs
2. **Configuration Analysis**: Evaluating anti-syn conformations using dynamic programming
3. **Linking Number**: Calculating optimal delta linking numbers
4. **Statistical Assessment**: Providing probability scores based on normal distribution

### Performance Optimizations

- **Parallel Processing**: Uses Rayon for multi-threaded computation
- **Memory Layout**: Optimized data structures for cache efficiency
- **Vectorization**: Loop unrolling and SIMD-friendly operations
- **Dynamic Programming**: Reduces exponential complexity to linear

## Dependencies

- [`clap`](https://crates.io/crates/clap): Command-line argument parsing
- [`anyhow`](https://crates.io/crates/anyhow): Error handling
- [`memmap2`](https://crates.io/crates/memmap2): Memory-mapped file I/O
- [`rand`](https://crates.io/crates/rand): Random number generation
- [`rayon`](https://crates.io/crates/rayon): Data parallelism

## Scientific Background

Z-DNA is a left-handed double helix form of DNA that can form under certain conditions, particularly in sequences with alternating purine-pyrimidine patterns. The formation of Z-DNA has biological significance in:

- Gene regulation
- Chromatin structure
- DNA-protein interactions
- Genomic instability

This tool helps researchers identify potential Z-DNA forming regions in genomic sequences, which can be important for understanding regulatory mechanisms and structural biology.

## Performance

The Rust implementation provides significant performance improvements over the original implementation:

- **Multi-threading**: Scales with available CPU cores
- **Memory Efficiency**: Reduced memory footprint through optimized data structures
- **Cache Optimization**: Improved memory access patterns
- **Vectorization**: Better utilization of modern CPU features

## File Formats

### Input

- **FASTA format**: Standard nucleotide sequence files
- **Supported bases**: A, T, G, C, N (case-insensitive)
- **Circular sequences**: Automatically handles sequence circularity

### Output

- **Z-SCORE files**: Tab-separated values with detailed analysis results
- **Progress reporting**: Real-time processing status updates

## Contributing

Contributions are welcome! Please ensure that:

1. Code follows Rust best practices
2. Performance optimizations maintain correctness
3. Tests are included for new features
4. Documentation is updated accordingly

## License

This project maintains the spirit of the original research while providing a modern, high-performance implementation.

## Citation

If you use Z-HUNT-3 in your research, please cite the original paper:

> Chou, P.J. and Ho, P.S. (1986) A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences. The EMBO Journal, 5(10), 2737-2744.

## Acknowledgments

- Original algorithm by Ping-jung Chou and Pui S. Ho
- Rust implementation optimizations and modernization
- Community contributions and feedback