
use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use std::path::Path;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io;


// Public writer function to write to gzip
pub fn writer(filename: &str) -> Box<dyn Write> {
    let path = Path::new(filename);
    let file = File::create(&path).unwrap();

    Box::new(BufWriter::with_capacity(
        128 * 1024,
        GzEncoder::new(file, Compression::default()),
    ))

}

// Reads in an input list of paths
pub fn read_input_list(filename: &str) -> Result<Vec<String>, io::Error>{

    let mut file = File::open(filename)?;
    let mut list = String::new();
    file.read_to_string(&mut list)?;

    let content: Vec<String> = list
        .split_whitespace()
        .map(|x| x.to_string())
        .collect();

    Ok(content)
}
