
#[cfg(test)]
mod tests {

    // use serial_test::serial;
    use crate::mol2::Mol2Reader;
    use crate::file_io::read_input_list;
    use crate::mol2utils;

    #[test]
    #[serial]
    fn read_mol2() {
        /*
        Tests whether a single mol2 can be read correctly
        */

        let filename = "data/test0000.mol2.gz";

        let mol2_reader = Mol2Reader::new(filename)
            .expect("Error: Failed Reading Test Data");

        let mut num_molecules = 0;
        for _mol in mol2_reader.into_iter() {
            num_molecules += 1;
        }

        assert!(num_molecules == 451);
    }

    #[test]
    #[serial]
    fn read_multi_mol2() {
        /*
        Tests whether a multiple mol2s can be read correctly
        */

        let input_files = vec![
            "data/test0000.mol2.gz",
            "data/test0001.mol2.gz",
            "data/test0002.mol2.gz",
            "data/test0003.mol2.gz",
            "data/test0004.mol2.gz",
            "data/test0005.mol2.gz"
        ];

        let mut num_molecules = 0;

        for filename in input_files.into_iter() {

            let mol2_reader = Mol2Reader::new(filename)
                .expect("Error: Failed Reading Test Data");

            for _mol in mol2_reader.into_iter() {
                num_molecules += 1;
            }

        };

        assert!(num_molecules == 6972);
    }

    #[test]
    #[serial]
    fn read_list() {
        /*
        Tests whether a list of input files can be read correctly
        */

        let input_list = "data/input_list.txt";
        let input_files = read_input_list(input_list)
            .expect("Error: Failed Reading Input List");

        let mut num_molecules = 0;
        for filename in input_files.iter() {

            let mol2_reader = Mol2Reader::new(filename)
                .expect("Error: Failed Reading Test Data");

            for _mol in mol2_reader.into_iter() {
                num_molecules += 1;
            }

        };

        assert!(num_molecules == 6972);
    }

    #[test]
    #[serial]
    fn run_grep_without_energy() {
        /*
        Tests whether a query can be read and processed without energy
        */

        let input_list = "data/input_list.txt";
        let input_files = read_input_list(input_list)
            .expect("Error: Failed Reading Input List");
        let output_filename = "test_grep_without_energy.mol2.gz";
        let query_filename = "data/zinc_list.txt";
        let tol = 1e-6;

        let num_passing = mol2utils::grep(
            input_files,
            query_filename,
            output_filename,
            tol
        ).unwrap();

        assert!(num_passing == 10);

    }

    #[test]
    #[serial]
    fn run_grep_with_energy() {

        let input_list = "data/input_list.txt";
        let input_files = read_input_list(input_list)
            .expect("Error: Failed Reading Input List");
        let output_filename = "test_grep_with_energy.mol2.gz";
        let query_filename = "data/zinc_list.tsv";
        let tol = 1e-6;

        let num_passing = mol2utils::grep(
            input_files,
            query_filename,
            output_filename,
            tol
        ).unwrap();

        assert!(num_passing == 8);
    }

    #[test]
    #[serial]
    fn run_split() {

        let input_list = "data/input_list.txt";
        let input_files = read_input_list(input_list)
            .expect("Error: Failed Reading Input List");
        let prefix = "split";
        let num_files = 10;

        let count_vec = mol2utils::split(
            input_files,
            prefix,
            num_files
        ).unwrap();

        let expected = vec![
            698, 698, 697, 697, 697,
            697, 697, 697, 697, 697
            ];

        
        assert!(expected.len() == count_vec.len());

        count_vec
            .into_iter()
            .enumerate()
            .for_each(|(idx, val)| {
                assert!(expected[idx] == val);
            });

    }

}
