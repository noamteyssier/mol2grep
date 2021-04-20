
use crate::mol2::Mol2;

use std::fs::File;
use std::io::Error;
use std::io::BufReader;
use std::io::prelude::*;

use std::collections::{HashMap, HashSet};


// Enumerate describing input query format
pub enum QueryFormat {
    WithScore(HashMap<Mol2, f64>),
    WithoutScore(HashSet<Mol2>)
}

// Struct describing file IO of input query
pub struct QueryReader {
    bufreader: BufReader<File>,
    line: String
}
impl QueryReader {

    // Inserts a molecule into a HashSet
    fn insert_to_set(&self, table: &mut HashSet<Mol2>, mol: Mol2) {
        table.insert(mol);
    }

    // Inserts a molecule into a HashMap
    fn insert_to_map(&self, table: &mut HashMap<Mol2, f64>, mol: Mol2) {
        let energy = mol.get_energy();
        table.insert(mol, energy);
    }

    // Creates a molecule with a given name
    fn mol_with_name(&self, name: &str) -> Mol2 {
        let mut mol = Mol2::new();
        mol.add_name(name.trim().to_owned());
        mol
    }

    // Creates a molecule with a given name and energy
    fn mol_with_name_and_energy(&self, name: &str, energy: &str) -> Mol2 {
        let mut mol = Mol2::new();
        mol.add_name(name.to_owned());
        mol.add_energy(
            energy.parse::<f64>()
            .expect("\n\nError: Malformed Energy Column...\n...Unable to be parsed into float\n\n"));
        mol
    }

    // Split a string on whitespace and return a vector of elements
    fn split_items(&self) -> Vec<&str> {
        self.line.trim()
            .split_whitespace()
            .collect()
    }

    // Read in a list of IDs without scores and construct a HashSet
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

    // Read in a list of IDS with scores and construct a HashMap
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

    // Step forward one line in file
    fn step(&mut self) -> Result<usize, Error> {
        /*
        Steps through the file one line at a time and returns false if EOF is found
        */
        self.line.clear();
        let eof = self.bufreader.read_line(&mut self.line)?;
        Ok(eof)
    }

    // Load in query input file with necessary format
    pub fn load_queries(&mut self) -> Result<QueryFormat, Error> {
        self.step()?;

        let items = self.split_items();

        match items.len() {
            1 => Ok(QueryFormat::WithoutScore(self.read_zinc_list())),
            2 => Ok(QueryFormat::WithScore(self.read_zinc_score_table())),
            _ => panic!("\n\nError: Malformed Query Input...\n..Found >2 columns but expecting 2\n\n")
        }
    }

    // Instantiate a new QueryReader
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
