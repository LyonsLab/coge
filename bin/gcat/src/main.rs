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
        writeln!(&mut stderr, "not enough arguments");
        process::exit(1);
    }

    let fasta_path = &args[1];
    let feat_type = &args[2];

    let fasta = read_fasta(fasta_path).unwrap();

    calc_percentages(&fasta);
}

fn calc_percentages(fasta: &str) {

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
