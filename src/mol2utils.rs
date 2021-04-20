
use std::thread;
use std::sync::mpsc;
use std::sync::mpsc::{Sender, Receiver};
use std::sync::{Arc, Mutex};

use std::collections::{HashMap, HashSet};
use std::io::Error;
use std::io::prelude::*;

use crate::mol2::{Mol2, Mol2Reader};
use crate::query::{QueryFormat, QueryReader};
use crate::file_io::writer;

use indicatif::ProgressIterator;
use rayon::prelude::*;

// Function to perform grep without checking for score matches
fn grep_with_set(
        mol2_reader: Mol2Reader,
        table: &HashSet<Mol2>,
        channel: &mut Sender<Mol2>) -> (u32, u32) {

    let mut num_molecules = 0;
    let mut num_passing = 0;

    mol2_reader
        .into_iter()
        .map(|x|{
            num_molecules += 1;
            x
        })
        .filter(|x|
            table.contains(x)
        )
        .for_each(|x|{
            num_passing += 1;
            channel.send(x).expect("Error: Broken Send Channel");
        });

    (num_molecules, num_passing)

}


// Function to perform grep while checking for score matches
fn grep_with_map(
        mol2_reader: Mol2Reader,
        table: &HashMap<Mol2, f64>,
        tol: f64,
        channel: &mut Sender<Mol2>) -> (u32, u32) {

    let mut num_molecules = 0;
    let mut num_passing = 0;

    mol2_reader
        .into_iter()
        .map(|x|{
            num_molecules += 1;
            x
        })
        .filter(|x|
            table.contains_key(x)
        )
        .filter(|x|
            (x.get_energy() - table.get(x).unwrap() <= tol)
        )
        .for_each(|x|{
            num_passing += 1;
            channel.send(x).expect("Error: Broken Send Channel");
        });

    (num_molecules, num_passing)

}


// implements grep subcommand
pub fn grep(
        input_files: Vec<String>,
        query_filename: &str,
        output_filename: &str,
        tol: f64) -> Result<u32, Error> {

    // Instantiate QueryReader and read file into table
    let mut qr = QueryReader::new(query_filename)?;
    let table = qr.load_queries()?;

    // Instantiate Writer
    let mut writer_file = writer(output_filename);

    // Instantiate Send/Receive Channels
    let (channel_send, channel_recv): (Sender<Mol2>, Receiver<Mol2>) = mpsc::channel();

    // Keep statistics on number of molecules processed
    let num_molecules = Arc::new(Mutex::new(0));
    let num_passing = Arc::new(Mutex::new(0));

    let num_molecules_fmt = num_molecules.clone();
    let num_passing_fmt = num_passing.clone();

    // places molecules into writer channel
    thread::spawn(move || {

        // iterate through input files in parallel
        input_files
            .into_iter()
            .progress()  // adds a progress bar on the file processing
            .par_bridge()
            .for_each_with(channel_send, |sender, x| {

                // instantiate a new mol2 reader
                let mol2_reader = Mol2Reader::new(&x).unwrap();

                // depending on the query input format
                let (nm, np) = match table {

                    // filter molecules without considering query score
                    QueryFormat::WithoutScore(ref t) => {
                        grep_with_set(mol2_reader, &t, sender)
                    },

                    // filter molecules considering query score
                    QueryFormat::WithScore(ref t) => {
                        grep_with_map(mol2_reader, &t, tol, sender)
                    }
                };

                *num_molecules.lock().unwrap() += nm;
                *num_passing.lock().unwrap() += np;

            });
    });

    // writes passing molecules to file
    for mol in channel_recv {
        writer_file
            .write_all(
                mol.get_lines().as_bytes()
            )
            .expect(
                "Error: Error writing to output file"
            )
    };

    println!(
        ">>> Number of Molecules Processed: {}",
        num_molecules_fmt.lock().unwrap()
    );

    println!(
        ">>> Number of Molecules Accepted: {}",
        num_passing_fmt.lock().unwrap()
    );

    let result = *num_passing_fmt
        .lock()
        .unwrap();

    Ok(result)
}

// implements split subcommand
pub fn split(
        input_files: Vec<String>,
        prefix: &str,
        num_files: usize) -> Result<Vec<u32>, Error> {

        // Instantiate Send/Receive Channels
        let (channel_send, channel_recv): (Sender<Mol2>, Receiver<Mol2>) = mpsc::channel();

        // places molecules into writer channel
        thread::spawn(move || {

            // iterate through input files in parallel
            input_files
                .into_iter()
                .progress()
                .par_bridge()
                .for_each_with(channel_send, |sender, x| {

                    // instantiate a new mol2 reader
                    let mol2_reader = Mol2Reader::new(&x).unwrap();

                    mol2_reader
                        .into_iter()
                        .for_each(|x|{
                            sender.send(x).expect("Error in sending through channel");
                        })

                });
        });

        let mut writer_vec: Vec<Box<dyn Write>> = (0..num_files)
            .into_iter()
            .map(|i| {
                writer(&format!("{}.{:04}.mol2.gz", prefix, i))
            })
            .collect();

        let mut count_vec = vec![0; num_files];


        let mut num_molecules = 0;

        for mol in channel_recv {

            let file_id = num_molecules % num_files;

            writer_vec[file_id]
                .write_all(mol.get_lines().as_bytes())
                .expect("Error in writing to output file");
            count_vec[file_id] += 1;

            num_molecules += 1;
        };

        println!("\nFile Totals:");
        (0..num_files)
            .into_iter()
            .for_each(|i| {
                println!("  {}.{:04}.mol2.gz:\t{}", prefix, i, count_vec[i])
            });

        Ok(count_vec)
}

pub fn table(
        input_files: Vec<String>,
        output_filename: &str,
        write_header: bool) -> Result<(), Error> {

    // Instantiate Send/Receive Channels
    let (channel_send, channel_recv): (Sender<Mol2>, Receiver<Mol2>) = mpsc::channel();

    // places molecules into writer channel
    thread::spawn(move || {

        // iterate through input files in parallel
        input_files
            .iter()
            .for_each(|x| {

                // instantiate a new mol2 reader
                let mol2_reader = Mol2Reader::new(&x).unwrap();

                mol2_reader
                    .into_iter()
                    .for_each(|x|{
                        channel_send.send(x).expect("Error in sending through channel");
                    })

            });
    });

    // Instantiate Writer
    let mut writer = writer(output_filename);

    // Writer a header if no_header flag isn't present
    if write_header {
        writer
            .write_all(
                &format!("ligand_id\tname\tenergy\n").into_bytes()
            )
            .expect("Error in writing to output file");
    }

    let mut ligand_id = 0;
    for mol in channel_recv {

        writer
            .write_all(
                &format!("{}\t{}\t{}\n", ligand_id, mol.get_name(), mol.get_energy()).into_bytes()
            )
            .expect("Error in writing to output file");

        ligand_id += 1;
    };

    println!("\n Total Poses: {}", ligand_id);
    println!(" Written to: {}", output_filename);

    Ok(())
}
