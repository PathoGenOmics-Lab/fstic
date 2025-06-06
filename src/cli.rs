use clap::{App, Arg, ArgMatches};
use rayon::ThreadPoolBuilder;
use std::path::PathBuf;

pub enum InputMode {
    Vcf,
    Table,
}

/// Builds the command-line interface using clap.
pub fn build_cli() -> App<'static, 'static> {
    App::new("fstic")
        .version("1.0.0")
        .about("Calculates pairwise genetic distances between samples using standard estimators.")
        .author("Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>")
        .arg(
            Arg::with_name("vcf")
                .long("vcf")
                .value_name("VCF_FILES")
                .help("One or more input VCF files.")
                .multiple(true)
                .required_unless_one(&["vcf_list", "table", "table_list"]),
        )
        .arg(
            Arg::with_name("vcf_list")
                .long("vcf-list")
                .value_name("VCF_LIST")
                .help("A file containing a list of input VCF files.")
                .required_unless_one(&["vcf", "table", "table_list"]),
        )
        .arg(
            Arg::with_name("table")
                .long("table")
                .value_name("TABLE_FILES")
                .help("One or more input table files (format detected by extension: .csv, .tsv, .tab).")
                .multiple(true)
                .required_unless_one(&["vcf", "vcf_list", "table_list"]),
        )
        .arg(
            Arg::with_name("table_list")
                .long("table-list")
                .value_name("TABLE_LIST")
                .help("A file containing a list of input table files.")
                .required_unless_one(&["vcf", "vcf_list", "table"]),
        )
        .arg(
            Arg::with_name("reference")
                .short("r")
                .long("reference")
                .value_name("FASTA_FILE")
                .help("Reference FASTA file (required for all inputs).")
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTPUT_FILE")
                .help("Output file name for the distance matrix.")
                .required(true),
        )
        .arg(
            Arg::with_name("formula")
                .short("f")
                .long("formula")
                .value_name("FORMULA")
                .help("The distance formula to use.")
                .possible_values(&["fst", "gst", "nei", "chord", "bray-curtis", "jost_d", "reynolds", "rogers"])
                .default_value("fst"),
        )
        .arg(
            Arg::with_name("normalize")
                .long("normalize")
                .help("Normalize by the number of loci (affects: fst, chord, bray-curtis, jost_d).")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("min_depth")
                .long("min-depth")
                .value_name("DP")
                .help("Minimum total read depth to keep a variant.")
                .default_value("30"),
        )
        .arg(
            Arg::with_name("min_af")
                .long("min-af")
                .value_name("FREQ")
                .help("Minimum alternate allele frequency to keep a variant.")
                .default_value("0.05"),
        )
        .arg(
            Arg::with_name("min_alt_reads")
                .long("min-alt-reads")
                .value_name("AD")
                .help("Minimum number of alternate allele reads to keep a variant.")
                .default_value("2"),
        )
        .arg(
            Arg::with_name("min_alt_rev_reads")
                .long("min-alt-rev-reads")
                .value_name("ADR")
                .help("Minimum number of alternate allele reverse reads to keep a variant.")
                .default_value("2"),
        )
        .arg(
            Arg::with_name("workers")
                .short("w")
                .long("workers")
                .value_name("NUM_WORKERS")
                .help("Number of worker threads (default: all available cores)."),
        )
}

pub fn get_input_files(matches: &ArgMatches) -> Result<(InputMode, Vec<PathBuf>), String> {
    if matches.is_present("vcf") || matches.is_present("vcf_list") {
        let files = get_files_from_args(matches, "vcf", "vcf_list")?;
        Ok((InputMode::Vcf, files))
    } else if matches.is_present("table") || matches.is_present("table_list") {
        let files = get_files_from_args(matches, "table", "table_list")?;
        Ok((InputMode::Table, files))
    } else {
        Err("No input files provided. Use --vcf, --vcf-list, --table, or --table-list.".to_string())
    }
}

fn get_files_from_args(matches: &ArgMatches, arg_direct: &str, arg_list: &str) -> Result<Vec<PathBuf>, String> {
    if let Some(list_file) = matches.value_of(arg_list) {
        super::io::read_file_list(list_file)
            .map_err(|e| format!("Error reading list file {}: {}", list_file, e))
    } else {
        Ok(matches
            .values_of(arg_direct)
            .unwrap()
            .map(PathBuf::from)
            .collect())
    }
}

pub fn configure_thread_pool(matches: &ArgMatches) {
    let num_workers = matches
        .value_of("workers")
        .and_then(|w| w.parse::<usize>().ok())
        .unwrap_or_else(num_cpus::get);

    ThreadPoolBuilder::new()
        .num_threads(num_workers)
        .build_global()
        .unwrap();
}
