
use std::fmt;
use std::hash::{Hash, Hasher};

use std::fs::File;
use std::io::Error;
use std::io::BufReader;
use std::io::prelude::*;

use flate2::read::MultiGzDecoder;
use regex::Regex;

// Struct representing molecular data from a mol2 formatted file
#[derive (Clone)]
pub struct Mol2 {
    name: String,
    energy: f64,
    lines: String
}
impl fmt::Debug for Mol2 {

    // Defines the debug formatting for a given Mol2
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Molecule")
         .field("\n  name", &self.name)
         .field("\n  energy", &self.energy)
         .finish()
    }

}
impl Hash for Mol2 {

    // Defines the hash of a Mol2 as the hash of its name
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }

}
impl PartialEq for Mol2 {

    // Defines equality between 2 Mol2 as based on Name
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}
impl Eq for Mol2 {}
impl Mol2 {

    // Instantiate a new Mol2
    pub fn new() -> Self {
        Mol2 {
            name: String::new(),
            energy: 100.0,
            lines: String::new()
        }
    }

    // Adds a name to current Mol2
    pub fn add_name(&mut self, name: String) {
        self.name = name;
    }

    // Adds an energy to current Mol2
    pub fn add_energy(&mut self, energy: f64) {
        self.energy = energy;
    }

    // Adds a raw text line to current Mol2
    pub fn add_line(&mut self, line: &str) {
        self.lines += line;
    }

    // Returns name from current Mol2
    pub fn get_name(&self) -> &str {
        return &self.name
    }

    // Returns energy from current Mol2
    pub fn get_energy(&self) -> f64 {
        self.energy
    }

    // Returns lines from current Mol2
    pub fn get_lines(&self) -> &str {
        return &self.lines
    }


}

// Struct describing file IO of a mol2 formatted file
pub struct Mol2Reader {
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

    // Instantiate a new Mol2Reader
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

    // Step forward one line in the file
    fn step(&mut self) -> bool {
        self.line.clear();
        let eof = self.reader.read_line(&mut self.line).unwrap();
        eof != 0
    }

    // Retrieve the next Mol2 in the file
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

            // Case where TRIPOS data is found
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
