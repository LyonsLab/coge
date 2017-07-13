use std::env;
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::prelude::*;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut stderr = std::io::stderr();

    if args.len() < 3 {
        writeln!(&mut stderr, "not enough arguments").unwrap();
        process::exit(1);
    }

    let fasta_path = &args[1];
//    let feat_type = &args[2];

    let fasta = read_fasta(fasta_path).unwrap();

    calc_percentages(&fasta);
}

fn calc_percentages(fasta: &str) {
    let bytes = fasta.as_bytes();
    let mut gc = 0;
    let mut at = 0;
    let mut nx = 0;
    let total = bytes.len();
    print!("{} ", total);
    for i in 0.. total {
        let c = bytes[i];
        if c == b'G' || c == b'g' || c == b'C' || c == b'c' {
            gc += 1;
        } else if c == b'A' || c == b'a' || c == b'T' || c == b't' {
            at += 1;
        } else if c == b'N' || c == b'n' || c == b'X' || c == b'x' {
            nx += 1;
        }
    }
    print!("{{\"gc\": {:.4}, \"at\": {:.4}, \"nx\": {:.4}}}", gc as f64 / total as f64 * 100.0, at  as f64/ total as f64 * 100.0, nx as f64 / total as f64 * 100.0);
}

// fn get_feature_indices() {

// }

fn read_fasta(fasta_path: &str) -> Result<String,io::Error> {
    let mut s = String::new();
    let f = File::open(fasta_path)?;
    let f = BufReader::new(f);
    for line in f.lines() {
        let l = line.unwrap();
        if !l.starts_with(">") {
            s += &l;
        }
    }
    Ok(s)
}
