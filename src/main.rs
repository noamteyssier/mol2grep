
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use io::Error;
use regex::Regex;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::io;
use std::fmt;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::path::Path;
use clap::{Arg, App, ArgMatches};


#[derive (Clone)]
struct Mol2 {
    name: String,
    energy: f64,
    lines: String
}
impl fmt::Debug for Mol2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Molecule")
         .field("\n  name", &self.name)
         .field("\n  energy", &self.energy)
         .finish()
    }
}
impl Hash for Mol2 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}
impl PartialEq for Mol2 {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}
impl Eq for Mol2 {}
impl Mol2 {
    pub fn new() -> Self {
        Mol2 {
            name: String::new(),
            energy: 100.0,
            lines: String::new()
        }
    }

    fn add_name(&mut self, name: String) {
        self.name = name;
    }

    fn add_energy(&mut self, energy: f64) {
        self.energy = energy;
    }

    fn add_line(&mut self, line: &str) {
        self.lines += line;
    }

}


struct Mol2Reader {
    reader: BufReader<MultiGzDecoder<File>>,
    line: String,
    regex_name: Regex,
    regex_energy: Regex,
    regex_tripos: Regex,
    regex_tripos_molecule: Regex
}
impl Iterator for Mol2Reader {
    type Item = Mol2;

    fn next(&mut self) -> Option<Mol2> {
        return self.get_mol2();
    }
}
impl Mol2Reader {
    pub fn new(filename: &str) -> Result<Self, Error> {
        let file = File::open(filename)?;
        let gzr = MultiGzDecoder::new(file);
        let reader = BufReader::new(gzr);
        let line = String::new();
        let regex_name = Regex::new(r"#+ +Name: +").unwrap();
        let regex_energy = Regex::new(r"#+ +Total Energy: +").unwrap();
        let regex_tripos = Regex::new(r"^@<TRIPOS>").unwrap();
        let regex_tripos_molecule = Regex::new(r"^@<TRIPOS>MOLECULE").unwrap();

        Ok(Mol2Reader {
            reader,
            line,
            regex_name,
            regex_energy,
            regex_tripos,
            regex_tripos_molecule
        })
    }

    fn step(&mut self) -> bool {
        self.line.clear();
        let eof = self.reader.read_line(&mut self.line).unwrap();
        eof != 0
    }

    fn get_mol2(&mut self) -> Option<Mol2> {
        let mut mol = Mol2::new();
        let mut tripos_state = 0;
        let mut tripos_counts = Vec::new();

        loop {

            if !self.step() {
                return None
            }

            // The beginning of a new molecule
            if self.regex_name.is_match(&self.line) {
                mol.add_name(
                    self.regex_name
                        .replace_all(&self.line, "")
                        .trim()
                        .to_string()
                );
            }

            // Adds the energy for a molecule
            else if self.regex_energy.is_match(&self.line) {
                mol.add_energy(
                    self.regex_energy
                        .replace_all(&self.line, "")
                        .trim()
                        .to_string()
                        .parse::<f64>()
                        .unwrap()
                )
            }

            else if self.regex_tripos.is_match(&self.line) {

                // case where TRIPOS Molecule Is found providing number of elements in following
                if self.regex_tripos_molecule.is_match(&self.line) {
                    mol.add_line(&self.line);

                    for _ in 0..2 {
                        if !self.step() {return None}
                        mol.add_line(&self.line);
                    }

                    tripos_counts = self.line
                        .trim()
                        .split_whitespace()
                        .map(|x| x.parse::<u8>().unwrap())
                        .collect();
                }

                else {
                    mol.add_line(&self.line);

                    for _ in 0..tripos_counts[tripos_state] {
                        if !self.step() {return None}
                        mol.add_line(&self.line);
                    }

                    tripos_state += 1;
                }

                self.line.clear();
            }

            if (tripos_state > 0) && (tripos_counts[tripos_state] == 0) {
                break;
            }

            mol.add_line(&self.line);
        }

        Some(mol)
    }
}


enum QueryFormat {
    WithScore(HashMap<Mol2, f64>),
    WithoutScore(HashSet<Mol2>)
}

struct QueryReader {
    bufreader: BufReader<File>,
    line: String
}
impl QueryReader {

    fn insert_to_set(&self, table: &mut HashSet<Mol2>, mol: Mol2) {
        table.insert(mol);
    }

    fn insert_to_map(&self, table: &mut HashMap<Mol2, f64>, mol: Mol2) {
        let energy = mol.energy;
        table.insert(mol, energy);
    }

    fn mol_with_name(&self, name: &str) -> Mol2 {
        let mut mol = Mol2::new();
        mol.add_name(name.trim().to_owned());
        mol
    }

    fn mol_with_name_and_energy(&self, name: &str, energy: &str) -> Mol2 {
        let mut mol = Mol2::new();
        mol.add_name(name.to_owned());
        mol.add_energy(
            energy.parse::<f64>()
            .expect("\n\nError: Malformed Energy Column...\n...Unable to be parsed into float\n\n"));
        mol
    }

    fn split_items(&self) -> Vec<&str> {
        self.line.trim()
            .split_whitespace()
            .collect()
    }

    fn read_zinc_list(&mut self) -> HashSet<Mol2> {
        let mut table = HashSet::new();

        let mol = self.mol_with_name(&self.line);
        self.insert_to_set(&mut table, mol);
        loop {
            if self.step().unwrap() == 0 {break;}

            let mol = self.mol_with_name(&self.line);
            self.insert_to_set(&mut table, mol);
        }

        table
    }

    fn read_zinc_score_table(&mut self) -> HashMap<Mol2, f64> {
        let mut table = HashMap::new();

        let items = self.split_items();

        let mol = self.mol_with_name_and_energy(items[0], items[1]);
        self.insert_to_map(&mut table, mol);

        loop {
            if self.step().unwrap() == 0 {break;}

            let items = self.split_items();
            let mol = self.mol_with_name_and_energy(items[0], items[1]);
            self.insert_to_map(&mut table, mol);
        }

        table
    }

    fn step(&mut self) -> Result<usize, Error> {
        /*
        Steps through the file one line at a time and returns false if EOF is found
        */
        self.line.clear();
        let eof = self.bufreader.read_line(&mut self.line)?;
        Ok(eof)
    }

    fn load_queries(&mut self) -> Result<QueryFormat, Error> {
        self.step()?;

        let items = self.split_items();

        match items.len() {
            1 => Ok(QueryFormat::WithoutScore(self.read_zinc_list())),
            2 => Ok(QueryFormat::WithScore(self.read_zinc_score_table())),
            _ => panic!("\n\nError: Malformed Query Input...\n..Found >2 columns but expecting 2\n\n")
        }
    }

    pub fn new(filename: &str) -> Result<Self, Error> {
        let file = File::open(filename)?;

        Ok(
            QueryReader {
                bufreader: BufReader::new(file),
                line: String::new()
            }
        )
    }
}


fn grep_with_set(mol2_reader: Mol2Reader, table: &HashSet<Mol2>, writer_file: &mut Box<dyn Write>) {

    mol2_reader
        .into_iter()
        .filter(|x| table.contains(x))
        .for_each(|x| writer_file.write_all(x.lines.as_bytes()).unwrap())

}

fn grep_with_map(mol2_reader: Mol2Reader, table: &HashMap<Mol2, f64>, tol: f64, writer_file: &mut Box<dyn Write>) {

    mol2_reader
        .into_iter()
        .filter(|x| table.contains_key(x))
        .filter(|x| (x.energy - table.get(x).unwrap() <= tol))
        .for_each(|x| writer_file.write_all(x.lines.as_bytes()).unwrap())

}

pub fn writer(filename: &str) -> Box<dyn Write> {
    let path = Path::new(filename);
    let file = File::create(&path).unwrap();

    Box::new(BufWriter::with_capacity(
        128 * 1024,
        GzEncoder::new(file, Compression::default()),
    ))

}


fn get_args() -> ArgMatches<'static> {
    let matches = App::new("mol2grep")
        .version("0.1")
        .author("Noam Teyssier")
        .about("Searches mol2 files and returns poses that match zincid and expected score")
        .arg(
            Arg::with_name("query")
                .short("q")
                .long("query")
                .value_name("ZINC-id,score.tsv")
                .help("Query table of ZINC-ids and scores to search for (tab separated, no header)")
                .takes_value(true)
                .required(true)
            )
        .arg(
            Arg::with_name("tolerance")
                .short("t")
                .long("tol")
                .help("Minimum tolerance to accept in float comparisons (default = 1e-6)")
                .takes_value(true)
                .required(false)
                .default_value("1e-6")
            )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("out")
                .help("mol2 formatted filename to write passing molecules to")
                .takes_value(true)
                .default_value("query_output.mol2.gz")
            )
        .arg(
            Arg::with_name("mol2")
                .short("i")
                .long("input")
                .value_name("*.mol2.gz")
                .help("mol2.gz formatted files to grep (can take multiple inputs)")
                .takes_value(true)
                .required(true)
                .min_values(1)
        )
        .get_matches();

    matches
}


fn main() -> Result<(), Error> {

    let matches = get_args();

    let query_filename = matches.value_of("query").unwrap();
    let output_filename = matches.value_of("output").unwrap();
    let input_files = matches.values_of("mol2").unwrap();
    let tol = matches
        .value_of("tolerance")
        .unwrap()
        .parse::<f64>()
        .unwrap();


    let mut qr = QueryReader::new(query_filename)?;
    let table = qr.load_queries()?;

    let mut writer_file = writer(output_filename);

    for filename in input_files {
        let mol2_reader = Mol2Reader::new(filename)?;

        match table {

            QueryFormat::WithoutScore(ref t) => {
                grep_with_set(mol2_reader, &t, &mut writer_file)
            },

            QueryFormat::WithScore(ref t) => {
                grep_with_map(mol2_reader, &t, tol, &mut writer_file)
            }
        }

    }

    Ok(())
}
