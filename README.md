# Mol2Grep

A simple tool for grepping a group of mol2 files with a specific list of queries (name, energy)

## Installation
```bash
git clone https://github.com/noamteyssier/mol2grep
cd mol2grep

cargo build --release
ln -s $(pwd)/target/release/mol2grep .

mol2grep --version
```

## Usage:
```bash
# example run 
mol2grep -i test0000.mol2.gz -q query_ids.tsv -o output.mol2.gz

# example run with multiple mol2 inputs
mol2grep -i test*.mol2.gz -q query_ids.tsv -o output.mol2.gz

# see options
mol2grep --help
```
