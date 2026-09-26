//! Reader for the `matrix/<case>/*.csv` goldens: a `row` column of quoted names, then
//! numeric columns; `NA` becomes `NaN`. Values are returned column-major so they drop
//! straight into the crate's matrix slices.

use std::env;
use std::path::{Path, PathBuf};

use super::CORPUS_ENV;

pub struct Matrix {
    pub file: String,
    pub row_names: Vec<String>,
    pub col_names: Vec<String>,
    pub nrow: usize,
    pub ncol: usize,
    /// Column-major `nrow x ncol`.
    pub data: Vec<f64>,
}

impl Matrix {
    pub fn read(dir: &Path, file: &str) -> Matrix {
        let path = dir.join(file);
        let mut reader = csv::Reader::from_path(&path)
            .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
        let headers: Vec<String> = reader
            .headers()
            .unwrap_or_else(|e| panic!("{file}: cannot read the header row: {e}"))
            .iter()
            .map(str::to_owned)
            .collect();
        assert_eq!(headers[0], "row", "{file}: first column must be `row`");
        let col_names = headers[1..].to_vec();
        let ncol = col_names.len();
        let mut row_names = Vec::new();
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for (i, record) in reader.records().enumerate() {
            let record = record.unwrap_or_else(|e| panic!("{file}: line {}: {e}", i + 2));
            row_names.push(record[0].to_owned());
            let vals: Vec<f64> = record
                .iter()
                .skip(1)
                .map(|cell| {
                    parse_cell(cell)
                        .unwrap_or_else(|| panic!("{file}: line {}: cannot parse {cell:?}", i + 2))
                })
                .collect();
            assert_eq!(vals.len(), ncol, "{file}: line {}: ragged row", i + 2);
            rows.push(vals);
        }
        let nrow = rows.len();
        let mut data = vec![0.0; nrow * ncol];
        for (i, row) in rows.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                data[j * nrow + i] = v;
            }
        }
        Matrix {
            file: file.to_owned(),
            row_names,
            col_names,
            nrow,
            ncol,
            data,
        }
    }

    pub fn col(&self, name: &str) -> &[f64] {
        let j = self
            .col_names
            .iter()
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("{}: no column {name:?}", self.file));
        &self.data[j * self.nrow..(j + 1) * self.nrow]
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.nrow + i]
    }
}

fn parse_cell(cell: &str) -> Option<f64> {
    match cell {
        "NA" | "NaN" => Some(f64::NAN),
        "Inf" => Some(f64::INFINITY),
        "-Inf" => Some(f64::NEG_INFINITY),
        "TRUE" => Some(1.0),
        "FALSE" => Some(0.0),
        _ => cell.parse().ok(),
    }
}

/// `matrix/<case>/` under the corpus root, or `None` with the skip printed.
pub fn matrix_dir(test: &str, case: &str) -> Option<PathBuf> {
    let root = env::var_os(CORPUS_ENV).filter(|value| !value.is_empty());
    let Some(root) = root.map(PathBuf::from) else {
        println!("{test}: skipped, {CORPUS_ENV} is unset or empty");
        return None;
    };
    if !root.is_dir() {
        println!(
            "{test}: skipped, {CORPUS_ENV}={} is not a directory",
            root.display()
        );
        return None;
    }
    let dir = root.join("matrix").join(case);
    assert!(
        dir.is_dir(),
        "{CORPUS_ENV}={} has no matrix/{case}/",
        root.display()
    );
    Some(dir)
}

/// The three matrix cases in the corpus.
pub const MATRIX_CASES: [&str; 3] = [
    "conditions_only",
    "conditions_numeric_covariate",
    "conditions_categorical_numeric_covariates",
];
