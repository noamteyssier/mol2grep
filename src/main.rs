
use clap::{Arg, App, ArgMatches, SubCommand, AppSettings};
use std::io::Error;

mod test;
mod mol2;
mod query;
mod mol2utils;
mod file_io;
use file_io::read_input_list;

// builds the global threadpool for rayon parallel processing
fn build_threadpool(num_threads: usize) {

    // Instantiate number of threads for rayon parallel processing
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Error : Failed to build thread pool");

}

// runs grep subcommand
fn subcommand_grep(matches: &ArgMatches) -> Result<(), Error> {

    // Assign Variables
    let input_mol2s = matches.values_of("mol2");
    let input_filelist = matches.value_of("input_files");

    let query_filename = matches.value_of("query").unwrap();
    let output_filename = matches.value_of("output").unwrap();
    let tol = matches.value_of("tolerance")
        .unwrap()
        .parse::<f64>()
        .expect("Malformed input: tolerance");

    let num_threads = matches.value_of("num_threads")
        .unwrap()
        .parse::<usize>()
        .expect("Malformed input: num_threads");

    // Instantiate Input File List
    let input_files: Vec<String>;
    match input_mol2s {

        // case where one or multiple mol2 are given at CLI
        Some(f) => {
            input_files = f.into_iter()
                .map(|x| x.to_string())
                .collect()
        },

        // case where a single input file containing mol2 paths is given at CLI
        None => {
            input_files = read_input_list(input_filelist.unwrap()).unwrap()
        }

    };

    build_threadpool(num_threads);

    mol2utils::grep(
        input_files,
        query_filename,
        output_filename,
        tol
    ).expect("Error: Failed to grep");

    Ok(())
}


// runs split subcommand
fn subcommand_split(matches: &ArgMatches) -> Result<(), Error> {

    // assign variables
    let input_mol2s = matches.values_of("mol2");
    let input_filelist = matches.value_of("input_files");
    let prefix = matches.value_of("prefix").unwrap();

    let num_files = matches.value_of("num_files")
        .unwrap()
        .parse::<usize>()
        .expect("Malformed input: num_threads");

    let num_threads = matches.value_of("num_threads")
        .unwrap()
        .parse::<usize>()
        .expect("Malformed input: num_threads");

    // Instantiate Input File List
    let input_files: Vec<String>;
    match input_mol2s {

        // case where one or multiple mol2 are given at CLI
        Some(f) => {
            input_files = f.into_iter()
                .map(|x| x.to_string())
                .collect()
        },

        // case where a single input file containing mol2 paths is given at CLI
        None => {
            input_files = read_input_list(input_filelist.unwrap()).unwrap()
        }

    };

    build_threadpool(num_threads);

    mol2utils::split(
        input_files,
        prefix,
        num_files
    ).expect("Error: Failed to split");

    Ok(())
}


// // runs table subcommand
fn subcommand_table(matches: &ArgMatches) -> Result<(), Error> {

    // assign variables
    let input_mol2s = matches.values_of("mol2");
    let input_filelist = matches.value_of("input_files");
    let output_filename = matches.value_of("output").unwrap();


    // Instantiate Input File List
    let input_files: Vec<String>;
    match input_mol2s {

        // case where one or multiple mol2 are given at CLI
        Some(f) => {
            input_files = f.into_iter()
                .map(|x| x.to_string())
                .collect()
        },

        // case where a single input file containing mol2 paths is given at CLI
        None => {
            input_files = read_input_list(input_filelist.unwrap()).unwrap()
        }

    };

    let add_header = !matches.is_present("no_header");

    mol2utils::table(
        input_files,
        output_filename,
        add_header
    )

}


// Receives arguments from CLI
fn build_cli() -> App<'static, 'static> {
    let app = App::new("mol2grep")
        .version("0.1")
        .author("Noam Teyssier")
        .subcommand(SubCommand::with_name("grep")
            .about("greps mol2 files and returns poses that match zincid and/or expected score")
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
                    .short("e")
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
                    .required_unless_one(&["input_files"])
            )
            .arg(
                Arg::with_name("input_files")
                .short("f")
                .long("files")
                .value_name("<files>.txt")
                .help("a list of filenames to process")
                .takes_value(true)
                .required(false)
            )
            .arg(
                Arg::with_name("num_threads")
                    .short("t")
                    .long("threads")
                    .help("Number of threads to use in parallel processing")
                    .takes_value(true)
                    .required(false)
                    .default_value("4")
                )
            .setting(AppSettings::ArgRequiredElseHelp)
        )
        .subcommand(SubCommand::with_name("split")
            .about("Splits a list of mol2 files into a given number of output files")
            .arg(
                Arg::with_name("mol2")
                    .short("i")
                    .long("input")
                    .value_name("*.mol2.gz")
                    .help("mol2.gz formatted files to grep (can take multiple inputs)")
                    .takes_value(true)
                    .required(true)
                    .min_values(1)
                    .required_unless_one(&["input_files"])
            )
            .arg(
                Arg::with_name("input_files")
                .short("f")
                .long("files")
                .value_name("<files>.txt")
                .help("a list of filenames to process")
                .takes_value(true)
                .required(false)
            )
            .arg(
                Arg::with_name("prefix")
                    .short("o")
                    .long("prefix")
                    .help("prefix of output files: <prefix>.file_id.mol2.gz")
                    .takes_value(true)
                    .default_value("split")
                )
            .arg(
                Arg::with_name("num_files")
                    .short("n")
                    .long("num_files")
                    .help("Number of files to split items into")
                    .takes_value(true)
                    .required(false)
                    .default_value("4")
                )
            .arg(
                Arg::with_name("num_threads")
                    .short("t")
                    .long("threads")
                    .help("Number of threads to use in parallel processing")
                    .takes_value(true)
                    .required(false)
                    .default_value("4")
                )
        )
        .subcommand(SubCommand::with_name("table")
            .about("convert a list of mol2 files into tab-separated table of names + scores")
            .arg(
                Arg::with_name("mol2")
                    .short("i")
                    .long("input")
                    .value_name("*.mol2.gz")
                    .help("mol2.gz formatted files to grep (can take multiple inputs)")
                    .takes_value(true)
                    .required(true)
                    .min_values(1)
                    .required_unless_one(&["input_files"])
            )
            .arg(
                Arg::with_name("input_files")
                .short("f")
                .long("files")
                .value_name("<files>.txt")
                .help("a list of filenames to process")
                .takes_value(true)
                .required(false)
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .help("output filename to write table to")
                    .takes_value(true)
                    .default_value("output.tab.gz")
                )
            .arg(
                Arg::with_name("no_header")
                    .short("n")
                    .long("no_header")
                    .help("do not include a header in output file")
                    .takes_value(false)
                )
        )
        .setting(AppSettings::SubcommandRequiredElseHelp);


    return app
}


fn main() {

    // Match Arguments
    let app = build_cli();
    let matches = app.get_matches();

    match matches.subcommand() {
        ("grep", grep_matches) => {
            subcommand_grep(grep_matches.unwrap())
                .expect("Error: Failed to grep")
        },
        ("split", split_matches) => {
            subcommand_split(split_matches.unwrap())
                .expect("Error: Failed to split")
        }
        ("table", table_matches) => {
            subcommand_table(table_matches.unwrap())
                .expect("Error: Failed to build table")
        }
        _ => unreachable!()
    };


}
