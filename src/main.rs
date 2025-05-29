//! Z-HUNT-3 computer program
//! 
//! Rust implementation of the Z-DNA formation prediction algorithm
//! Originally written by Ping-jung Chou, under the instruction of Pui S. Ho
//! Based on the paper "A computer aided thermodynamic approach for predicting 
//! the formation of Z-DNA in naturally occurring sequences"
//! The EMBO Journal, Vol.5, No.10, pp2737-2744, 1986

use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::{SystemTime, UNIX_EPOCH};
use clap::Parser;
use anyhow::{Result, Context};
use rayon::prelude::*;
use std::sync::Arc;

/// Z-HUNT: Predicts Z-DNA formation in DNA sequences
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Window size (dinucleotides)
    window_size: usize,
    
    /// Minimum size
    min_size: usize,
    
    /// Maximum size
    max_size: usize,
    
    /// Input sequence file
    filename: String,
    
    /// Number of threads to use (0 = auto-detect)
    #[arg(short, long, default_value = "0")]
    threads: usize,
    
    /// Use sequential processing instead of parallel
    #[arg(short, long)]
    sequential: bool,
}

// Constants
const K_RT: f64 = -0.2521201;  // -1100/4363
const SIGMA: f64 = 16.94800353; // 10/RT
const EXP_LIMIT: f64 = -600.0;
const RT: f64 = 0.59004; // 0.00198*298
const A: f64 = 0.357; // 2 * (1/10.5 + 1/12)
const B: f64 = 0.4;

/// Format probability with leading zero in exponent (e.g., 3.026373e-02 instead of 3.026373e-2)
/// Special case: zero exponent should be formatted as e0 (not e+00)
fn format_probability(prob: f64) -> String {
    let formatted = format!("{:.6e}", prob);
    
    // Check if we have a single-digit exponent and add leading zero
    if let Some(e_pos) = formatted.find('e') {
        let (mantissa, exponent_part) = formatted.split_at(e_pos);
        let exponent_str = &exponent_part[1..]; // Skip the 'e'
        
        if exponent_str.len() == 2 && (exponent_str.starts_with('+') || exponent_str.starts_with('-')) {
            let sign = &exponent_str[0..1];
            let digit = &exponent_str[1..2];
            
            // Special case: zero exponent should be formatted as e0
            if digit == "0" {
                format!("{}e0", mantissa)
            } else {
                // Single digit exponent, add leading zero
                format!("{}e{}0{}", mantissa, sign, digit)
            }
        } else {
            formatted
        }
    } else {
        formatted
    }
}

// Delta BZ Energy of Dinucleotide
const DBZED: [[f64; 17]; 4] = [
    // AS-AS
    [4.40, 6.20, 3.40, 5.20, 2.50, 4.40, 1.40, 3.30, 3.30, 5.20, 2.40, 4.20, 1.40, 3.40, 0.66, 2.40, 4.26],
    // SA-SA  
    [4.40, 2.50, 3.30, 1.40, 6.20, 4.40, 5.20, 3.40, 3.40, 1.40, 2.40, 0.66, 5.20, 3.30, 4.20, 2.40, 4.26],
    // AS-SA
    [6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00, 4.26],
    // SA-AS
    [6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00, 4.26],
];

#[derive(Clone)]
struct ZHunt {
    exp_dbzed: [[f64; 17]; 4],
    bzindex: Vec<usize>,
    bzenergy: Vec<f64>,
    best_bzenergy: Vec<f64>,
    bztwist: Vec<f64>,
    logcoef: Vec<f64>,
    antisyn: Vec<char>,
    best_antisyn: Vec<char>,
    best_esum: f64,
    deltatwist: f64,
    terms: usize,
}

#[derive(Clone)]
struct PositionResult {
    position: usize,
    end_position: usize,
    length: usize,
    best_dl: f64,
    slope: f64,
    probability: f64,
    sequence_segment: String,
    antisyn_str: String,
}

impl ZHunt {
    fn new(dinucleotides: usize) -> Self {
        let nucleotides = 2 * dinucleotides;
        
        // Calculate exp(-dbzed/rt) - pre-computed for better cache performance
        let mut exp_dbzed = [[0.0; 17]; 4];
        for i in 0..4 {
            for j in 0..17 {
                exp_dbzed[i][j] = (-DBZED[i][j] / RT).exp();
            }
        }
        
        // Calculate bztwist values with better memory layout
        let mut bztwist = Vec::with_capacity(dinucleotides);
        let mut ab = B + B;
        for _ in 0..dinucleotides {
            ab += A;
            bztwist.push(ab);
        }
        
        Self {
            exp_dbzed,
            bzindex: vec![0; dinucleotides],
            bzenergy: vec![0.0; dinucleotides],
            best_bzenergy: vec![0.0; dinucleotides],
            bztwist,
            logcoef: vec![0.0; dinucleotides],
            antisyn: vec![' '; nucleotides + 1],
            best_antisyn: vec![' '; nucleotides + 1],
            best_esum: 0.0,
            deltatwist: 0.0,
            terms: 0,
        }
    }
    
    #[inline]
    fn linear_search<F>(&self, x1: f64, x2: f64, tolerance: f64, func: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let f = func(x1);
        let fmid = func(x2);
        
        if f * fmid >= 0.0 {
            return x2;
        }
        
        let (mut x, mut dx) = if f < 0.0 {
            (x1, x2 - x1)
        } else {
            (x2, x1 - x2)
        };
        
        // Limit iterations to prevent infinite loops
        for _ in 0..100 {
            dx *= 0.5;
            let xmid = x + dx;
            let fmid = func(xmid);
            if fmid <= 0.0 {
                x = xmid;
            }
            if dx.abs() <= tolerance {
                break;
            }
        }
        x
    }
    
    #[inline]
    fn delta_linking(&self, dl: f64) -> f64 {
        let mut expmini = 0.0;
        
        // Find minimum exponent with better vectorization
        let bztwist = &self.bztwist[..self.terms];
        let logcoef = &self.logcoef[..self.terms];
        
        for i in 0..self.terms {
            let z = dl - bztwist[i];
            let exp_val = logcoef[i] + K_RT * z * z;
            if exp_val < expmini {
                expmini = exp_val;
            }
        }
        
        expmini = if expmini < EXP_LIMIT {
            EXP_LIMIT - expmini
        } else {
            0.0
        };
        
        let mut sump = 0.0;
        let mut sumq = 0.0;
        
        // Unroll loop for better performance
        let mut i = 0;
        while i + 3 < self.terms {
            for j in 0..4 {
                let z = dl - bztwist[i + j];
                let exp_val = (logcoef[i + j] + K_RT * z * z + expmini).exp();
                sumq += exp_val;
                sump += bztwist[i + j] * exp_val;
            }
            i += 4;
        }
        
        // Handle remaining elements
        while i < self.terms {
            let z = dl - bztwist[i];
            let exp_val = (logcoef[i] + K_RT * z * z + expmini).exp();
            sumq += exp_val;
            sump += bztwist[i] * exp_val;
            i += 1;
        }
        
        sumq += (K_RT * dl * dl + SIGMA + expmini).exp();
        self.deltatwist - sump / sumq
    }
    
    #[inline]
    fn delta_linking_slope(&self, dl: f64) -> f64 {
        let mut expmini = 0.0;
        
        // Find minimum exponent with better vectorization
        let bztwist = &self.bztwist[..self.terms];
        let logcoef = &self.logcoef[..self.terms];
        
        for i in 0..self.terms {
            let z = dl - bztwist[i];
            let exp_val = logcoef[i] + K_RT * z * z;
            if exp_val < expmini {
                expmini = exp_val;
            }
        }
        
        expmini = if expmini < EXP_LIMIT {
            EXP_LIMIT - expmini
        } else {
            0.0
        };
        
        let mut sump = 0.0;
        let mut sump1 = 0.0;
        let mut sumq = 0.0;
        let mut sumq1 = 0.0;
        let x = 2.0 * K_RT;
        
        // Unroll loop for better performance
        let mut i = 0;
        while i + 3 < self.terms {
            for j in 0..4 {
                let z = dl - bztwist[i + j];
                let y = (logcoef[i + j] + K_RT * z * z + expmini).exp();
                sumq += y;
                sump += bztwist[i + j] * y;
                let y_scaled = y * z * x;
                sumq1 += y_scaled;
                sump1 += bztwist[i + j] * y_scaled;
            }
            i += 4;
        }
        
        // Handle remaining elements
        while i < self.terms {
            let z = dl - bztwist[i];
            let y = (logcoef[i] + K_RT * z * z + expmini).exp();
            sumq += y;
            sump += bztwist[i] * y;
            let y_scaled = y * z * x;
            sumq1 += y_scaled;
            sump1 += bztwist[i] * y_scaled;
            i += 1;
        }
        
        let y = (K_RT * dl * dl + SIGMA + expmini).exp();
        sumq += y;
        sumq1 += x * dl * y;
        
        (sump1 - sump * sumq1 / sumq) / sumq
    }
    
    #[inline]
    fn assign_bzenergy_index(&mut self, nucleotides: usize, sequence: &[char]) {
        let mut j = 0;
        let mut i = 0;
        
        while i < nucleotides {
            let c1 = sequence[i];
            let c2 = sequence[i + 1];
            
            // Optimized lookup using direct byte comparison
            let idx = match (c1 as u8, c2 as u8) {
                (b'a', b'a') => 0,
                (b'a', b't') => 1,
                (b'a', b'g') => 2,
                (b'a', b'c') => 3,
                (b't', b'a') => 4,
                (b't', b't') => 5,
                (b't', b'g') => 6,
                (b't', b'c') => 7,
                (b'g', b'a') => 8,
                (b'g', b't') => 9,
                (b'g', b'g') => 10,
                (b'g', b'c') => 11,
                (b'c', b'a') => 12,
                (b'c', b't') => 13,
                (b'c', b'g') => 14,
                (b'c', b'c') => 15,
                _ => 16, // Any combination with 'n' or unknown
            };
            
            self.bzindex[j] = idx;
            j += 1;
            i += 2;
        }
    }
    
    
    fn anti_syn_energy(&mut self, _din: usize, dinucleotides: usize, _esum: f64) {
        // Use dynamic programming instead of exponential recursion
        // This reduces complexity from O(2^n) to O(n)
        
        // For each position and configuration (AS=0, SA=1), store the minimum energy
        let mut dp_energy = vec![vec![f64::INFINITY; 2]; dinucleotides + 1];
        let mut dp_bzenergy = vec![vec![Vec::<f64>::new(), Vec::<f64>::new()]; dinucleotides + 1];
        let mut dp_antisyn = vec![vec![Vec::<char>::new(), Vec::<char>::new()]; dinucleotides + 1];
        
        // Base case: empty sequence
        dp_energy[0][0] = 0.0;
        dp_energy[0][1] = 0.0;
        
        // Fill DP table
        for din in 0..dinucleotides {
            let bzindex_din = self.bzindex[din];
            
            // Try both previous configurations
            for prev_config in 0..2 {
                if dp_energy[din][prev_config] == f64::INFINITY {
                    continue;
                }
                
                // Try AS configuration for current position
                let i = if din == 0 {
                    0
                } else if prev_config == 0 { // previous was AS (ends with 'S')
                    0
                } else { // previous was SA (ends with 'A')
                    3
                };
                
                let e_as = DBZED[i][bzindex_din];
                let new_energy_as = dp_energy[din][prev_config] + e_as;
                
                if new_energy_as < dp_energy[din + 1][0] {
                    dp_energy[din + 1][0] = new_energy_as;
                    
                    dp_bzenergy[din + 1][0] = dp_bzenergy[din][prev_config].clone();
                    dp_bzenergy[din + 1][0].push(self.exp_dbzed[i][bzindex_din]);
                    
                    dp_antisyn[din + 1][0] = dp_antisyn[din][prev_config].clone();
                    dp_antisyn[din + 1][0].push('A');
                    dp_antisyn[din + 1][0].push('S');
                }
                
                // Try SA configuration for current position
                let i = if din == 0 {
                    1
                } else if prev_config == 1 { // previous was SA (ends with 'A')
                    1
                } else { // previous was AS (ends with 'S')
                    2
                };
                
                let e_sa = DBZED[i][bzindex_din];
                let new_energy_sa = dp_energy[din][prev_config] + e_sa;
                
                if new_energy_sa < dp_energy[din + 1][1] {
                    dp_energy[din + 1][1] = new_energy_sa;
                    
                    dp_bzenergy[din + 1][1] = dp_bzenergy[din][prev_config].clone();
                    dp_bzenergy[din + 1][1].push(self.exp_dbzed[i][bzindex_din]);
                    
                    dp_antisyn[din + 1][1] = dp_antisyn[din][prev_config].clone();
                    dp_antisyn[din + 1][1].push('S');
                    dp_antisyn[din + 1][1].push('A');
                }
            }
        }
        
        // Find the best configuration
        let best_config = if dp_energy[dinucleotides][0] < dp_energy[dinucleotides][1] { 0 } else { 1 };
        
        // Update the best solution
        self.best_esum = dp_energy[dinucleotides][best_config];
        self.best_bzenergy[..dinucleotides].copy_from_slice(&dp_bzenergy[dinucleotides][best_config]);
        
        // Copy antisyn configuration
        let best_antisyn_vec = &dp_antisyn[dinucleotides][best_config];
        for (i, &ch) in best_antisyn_vec.iter().enumerate() {
            if i < self.best_antisyn.len() {
                self.best_antisyn[i] = ch;
            }
        }
        
        // Null terminate
        if best_antisyn_vec.len() < self.best_antisyn.len() {
            self.best_antisyn[best_antisyn_vec.len()] = '\0';
        }
    }
    
    fn find_delta_linking(&mut self, dinucleotides: usize) -> f64 {
        // Initialize bzenergy with fill
        self.bzenergy[..dinucleotides].fill(1.0);
        
        // Calculate logcoef with better memory access pattern
        for i in 0..dinucleotides {
            let mut sum = 0.0;
            let remaining = dinucleotides - i;
            
            // Vectorized multiplication and sum
            for j in 0..remaining {
                self.bzenergy[j] *= self.best_bzenergy[i + j];
                sum += self.bzenergy[j];
            }
            self.logcoef[i] = sum.ln();
        }
        
        self.terms = dinucleotides;
        self.linear_search(10.0, 50.0, 0.001, |dl| self.delta_linking(dl))
    }
    
    #[inline]
    fn assign_probability(&self, dl: f64) -> f64 {
        const AVERAGE: f64 = 29.6537135;
        const STDV: f64 = 2.71997;
        const SQRT2_INV: f64 = 0.70710678118654752440; // 1/sqrt(2)
        const SQRTPI_INV: f64 = 0.564189583546; // 1/sqrt(pi)
        
        let z = (dl - AVERAGE).abs() / STDV;
        let mut x = z * SQRT2_INV;
        let y = SQRTPI_INV * (-x * x).exp();
        let z_sq = z * z;
        let mut k = 1.0;
        let mut sum = 0.0;
        
        // Unroll first few iterations for better performance
        for _ in 0..8 {
            let old_sum = sum;
            sum += x;
            k += 2.0;
            x *= z_sq / k;
            if sum + x <= old_sum {
                break;
            }
        }
        
        // Continue with regular loop if needed
        loop {
            let old_sum = sum;
            sum += x;
            k += 2.0;
            x *= z_sq / k;
            if sum + x <= old_sum || k > 100.0 {
                break;
            }
        }
        
        let tail_prob = 0.5 - y * sum;
        if dl > AVERAGE {
            tail_prob
        } else {
            1.0 / tail_prob
        }
    }

    fn process_position(
        &mut self,
        position: usize,
        sequence: &[char],
        nucleotides: usize,
        from_din: usize,
        to_din: usize,
        a_half: f64,
        init_esum: f64,
    ) -> PositionResult {
        const PI_DEG: f64 = 57.29577951; // 180/pi
        
        self.assign_bzenergy_index(nucleotides, &sequence[position..]);
        let mut best_dl = 50.0;
        let mut best_antisyn = vec![' '; nucleotides + 1];

        for din in from_din..=to_din {
            self.best_esum = init_esum;
            self.deltatwist = a_half * din as f64;
            self.antisyn[2 * din] = '\0';
            self.anti_syn_energy(0, din, 0.0);
            let dl = self.find_delta_linking(din);

            if dl < best_dl {
                best_dl = dl;
                best_antisyn.copy_from_slice(&self.best_antisyn);
            }
        }

        let antisyn_str: String = best_antisyn.iter().take_while(|&&c| c != '\0').collect();
        let j = antisyn_str.len();
        let slope = self.delta_linking_slope(best_dl).atan() * PI_DEG;
        let probability = self.assign_probability(best_dl);

        let sequence_segment: String = sequence[position..position + j].iter().collect();

        PositionResult {
            position: position + 1, // 1-based indexing for output
            end_position: position + 1 + j,
            length: j,
            best_dl,
            slope,
            probability,
            sequence_segment,
            antisyn_str,
        }
    }
}

fn input_sequence(filename: &str, nucleotides: usize) -> Result<Vec<char>> {
    use std::fs;
    
    println!("Inputting sequence from {}", filename);
    
    // Read entire file at once for better performance
    let contents = fs::read_to_string(filename)
        .with_context(|| format!("Failed to read file: {}", filename))?;
    
    // Pre-allocate with better estimate
    let mut all_chars = Vec::with_capacity(contents.len());
    
    // Process bytes directly for better performance
    for &byte in contents.as_bytes() {
        match byte {
            b'a' | b'A' => all_chars.push('a'),
            b't' | b'T' => all_chars.push('t'),
            b'g' | b'G' => all_chars.push('g'),
            b'c' | b'C' => all_chars.push('c'),
            b'n' | b'N' => all_chars.push('n'),
            _ => {} // Skip non-nucleotide characters
        }
    }
    
    let original_len = all_chars.len();
    println!("Loaded {} nucleotides", original_len);
    
    // Add circular nucleotides
    all_chars.reserve(nucleotides);
    for j in 0..nucleotides {
        if j < original_len {
            all_chars.push(all_chars[j]);
        }
    }
    
    Ok(all_chars)
}

fn calculate_zscore_parallel(
    zhunt: &mut ZHunt,
    a: f64,
    max_dinucleotides: usize,
    min: usize,
    max: usize,
    filename: &str,
) -> Result<()> {
    let from_din = min;
    let to_din = max.min(max_dinucleotides);
    let nucleotides = 2 * to_din;
    
    println!("Calculating zscore with multithreading");
    
    let sequence = Arc::new(input_sequence(filename, nucleotides)?);
    let seq_length = sequence.len() - nucleotides;
    
    let output_filename = format!("{}.Z-SCORE", filename);
    let output_file = File::create(&output_filename)
        .with_context(|| format!("Failed to create output file: {}", output_filename))?;
    let mut writer = BufWriter::with_capacity(64 * 1024, output_file); // Larger buffer
    
    writeln!(writer, "{} {} {} {}", filename, seq_length, from_din, to_din)?;
    
    let a_half = a / 2.0;
    let init_esum = 10.0 * to_din as f64;
    
    let start_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_micros();
    
    // Optimize chunk size for better cache locality and load balancing
    let num_threads = rayon::current_num_threads();
    let chunk_size = (seq_length / (num_threads * 8)).max(1000).min(5000);
    
    for chunk_start in (0..seq_length).step_by(chunk_size) {
        let chunk_end = (chunk_start + chunk_size).min(seq_length);
        let sequence_ref = Arc::clone(&sequence);
        
        let mut chunk_results: Vec<PositionResult> = (chunk_start..chunk_end)
            .into_par_iter()
            .map(|i| {
                let mut local_zhunt = zhunt.clone();
                local_zhunt.process_position(
                    i,
                    &sequence_ref,
                    nucleotides,
                    from_din,
                    to_din,
                    a_half,
                    init_esum,
                )
            })
            .collect();
        
        // Sort chunk results by position
        chunk_results.sort_unstable_by_key(|r| r.position);
        
        // Write chunk results immediately to reduce memory usage
        for result in chunk_results {
            write!(writer, "{} {} {} {:7.3} {:7.3} {} ",
                   result.position, result.end_position, result.length,
                   result.best_dl, result.slope, format_probability(result.probability))?;
            
            write!(writer, "{}   {}", result.sequence_segment, result.antisyn_str)?;
            writeln!(writer)?;
        }
        
        // Progress indicator - less frequent updates
        if chunk_end % 20000 == 0 || chunk_end == seq_length {
            println!("Processed {}/{} positions", chunk_end, seq_length);
        }
    }
    
    let end_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_micros();
    
    println!("Run time: {} micros", end_time - start_time);
    Ok(())
}

fn calculate_zscore_sequential(
    zhunt: &mut ZHunt,
    a: f64,
    max_dinucleotides: usize,
    min: usize,
    max: usize,
    filename: &str,
) -> Result<()> {
    let from_din = min;
    let to_din = max.min(max_dinucleotides);
    let nucleotides = 2 * to_din;
    
    println!("Calculating zscore sequentially");
    
    let sequence = input_sequence(filename, nucleotides)?;
    let seq_length = sequence.len() - nucleotides;
    
    let output_filename = format!("{}.Z-SCORE", filename);
    let output_file = File::create(&output_filename)
        .with_context(|| format!("Failed to create output file: {}", output_filename))?;
    let mut writer = BufWriter::new(output_file);
    
    writeln!(writer, "{} {} {} {}", filename, seq_length, from_din, to_din)?;
    
    let a_half = a / 2.0;
    let init_esum = 10.0 * to_din as f64;
    
    let start_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();
    
    // Process all positions sequentially
    for i in 0..seq_length {
        let result = zhunt.process_position(
            i,
            &sequence,
            nucleotides,
            from_din,
            to_din,
            a_half,
            init_esum,
        );
        
        write!(writer, "{} {} {} {:7.3} {:7.3} {} ",
               result.position, result.end_position, result.length,
               result.best_dl, result.slope, format_probability(result.probability))?;
        
        write!(writer, "{}   {}", result.sequence_segment, result.antisyn_str)?;
        writeln!(writer)?;
    }
    
    let end_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();
    
    println!("Run time: {} sec", end_time - start_time);
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    println!("dinucleotides {}", args.window_size);
    println!("min/max {} {}", args.min_size, args.max_size);
    println!("operating on {}", args.filename);
    
    // Configure thread pool
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .with_context(|| format!("Failed to configure {} threads", args.threads))?;
        println!("Using {} threads", args.threads);
    } else {
        let num_cpus = rayon::current_num_threads();
        println!("Using {} threads (auto-detected)", num_cpus);
    }
    
    let mut zhunt = ZHunt::new(args.window_size);
    
    if args.sequential {
        calculate_zscore_sequential(
            &mut zhunt,
            A,
            args.window_size,
            args.min_size,
            args.max_size,
            &args.filename,
        )?;
    } else {
        calculate_zscore_parallel(
            &mut zhunt,
            A,
            args.window_size,
            args.min_size,
            args.max_size,
            &args.filename,
        )?;
    }
    
    Ok(())
}
