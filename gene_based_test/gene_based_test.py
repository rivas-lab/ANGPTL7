from __future__ import division
import argparse


def is_pos_def_and_full_rank(X):

    """ 
    Ensures a matrix is positive definite and full rank.
  
    Keep diagonals and multiples every other cell by .99.
  
    Parameters: 
    X: Matrix to verify.
  
    Returns: 
    X: Verified (and, if applicable, adjusted) matrix.
  
    """
    # Check that it's not only pos def but also above a certain threshold
    i = 0
    while not (np.all(np.linalg.eigvals(X) > 1e-15)):
        X = 0.99 * X + 0.01 * np.diag(np.diag(X))
        i += 1
        if i > 5:
            return np.nan
    return X


def safe_inv(X, matrix_name, block):

    """ 
    Safely inverts a matrix, or returns NaN.
  
    Parameters: 
    X: Matrix to invert.
    matrix_name: One of "U"/"v_beta" - used to print messages when inversion fails.
    block: Name of the aggregation block (gene). 
        Used to print messages when inversion fails.
  
    Returns: 
    X_inv: Inverse of X.
  
    """

    try:
        X_inv = np.linalg.inv(X)
    except LinAlgError as err:
        print(
            "Could not invert " + matrix_name + " for gene " + block + "."
        )
        return np.nan
    return X_inv


def farebrother(quad_T, d, fb):

    """ 
    Farebrother method from CompQuadForm.
  
    Parameters: 

    quad_T: Value point at which distribution function is to be evaluated.
    d: Distinct non-zero characteristic root(s) of A*Sigma. 
    fb: Farebrother R method (rpy2 object).
  
    Returns: 
    p_value: Farebrother p-value.
  
    """

    res = fb(quad_T, d)
    return np.asarray(res)[0]


def initialize_r_objects():

    """ 
    Initializes Farebrother R method as rpy2 object.
  
    Returns: 
    fb: Farebrother R method (rpy2 object).
  
    """

    robjects.r(
        """
    require(MASS)
    require(CompQuadForm)
    farebrother.method <- function(quadT, d, h = rep(1, length(d)), delta = rep(0, length(d)), maxiter = 100000, epsilon = 10^-16, type = 1) {
        return(farebrother(quadT, d, h, delta, maxit = as.numeric(maxiter), eps = as.numeric(epsilon), mode = as.numeric(type))$Qq)
    }
    """
    )
    fb = robjects.globalenv["farebrother.method"]
    fb = robjects.r["farebrother.method"]
    return fb


def return_BF_pval(U, beta, v_beta, fb, gene):

    """ 
    Computes a p-value from the quadratic form that is subsumed by the Bayes Factor.
  
    Parameters: 

    U: Kronecker product of the three matrices (S*M*K x S*M*K)
        dictating correlation structures; no missing data.
    beta: Effect size vector without missing data.
    v_beta: Diagonal matrix of variances of effect sizes without missing data.
    fb: Farebrother R method (rpy2 object).
    gene: Name of the gene.
  
    Returns: 
    p_value: Farebrother p-value.
  
    """
    v_beta = is_pos_def_and_full_rank(v_beta)
    if np.any(np.isnan(v_beta)):
        return [np.nan]
    U = is_pos_def_and_full_rank(U)
    if np.any(np.isnan(U)):
        return [np.nan]
    v_beta_inv = safe_inv(v_beta, "v_beta", gene)
    n = beta.shape[0]
    A = v_beta + U
    A = is_pos_def_and_full_rank(A)
    if np.any(np.isnan(A)):
        return [np.nan]
    A_inv = np.linalg.inv(A)
    quad_T = np.asmatrix(beta.T) * np.asmatrix((v_beta_inv - A_inv)) * np.asmatrix(beta)
    B = is_pos_def_and_full_rank(npm.eye(n) - np.asmatrix(A_inv) * np.asmatrix(v_beta))
    if np.any(np.isnan(B)):
        return [np.nan]
    d = np.linalg.eig(B)[0]
    d = [i for i in d if i > 0.01]
    p_value = farebrother(quad_T, d, fb)
    p_value = max(0, min(1, p_value))
    return [p_value]


def delete_rows_and_columns(X, indices_to_remove):

    """ 
    Helper function to delete rows and columns from a matrix.
  
    Parameters: 
    X: Matrix that needs adjustment.
    indices_to_remove: Rows and columns to be deleted.
  
    Returns: 
    X: Smaller matrix that has no missing data.
  
    """

    X = np.delete(X, indices_to_remove, axis=0)
    X = np.delete(X, indices_to_remove, axis=1)
    return X


def adjust_for_missingness(U, omega, beta, se, beta_list):

    """ 
    Deletes rows and columns where we do not have effect sizes/standard errors.

    Calls method delete_rows_and_columns, a helper function that calls the numpy
        command.
  
    Parameters: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K) 
        dictating correlation structures; may relate to missing data.
    omega: (S*M*K x S*M*K) matrix that contains correlation of errors
        across variants, studies, and phenotypes. 
    beta: Vector of effect sizes within the unit of aggregation;
        may contain missing data.
    se: Vector of standard errors within the unit of aggregation;
        may contain missing data.
    beta_list: List of effect sizes within the unit of aggregation;
        may contain missing data.
  
    Returns: 
    U: Potentially smaller U matrix not associated with missing data.
    omega: Potentially smaller omega matrix not associated with missing data.
    beta: Potentially smaller beta vector without missing data.
    se: Potentially smaller SE vector without missing data.
  
    """

    indices_to_remove = np.argwhere(np.isnan(beta_list))
    U = delete_rows_and_columns(U, indices_to_remove)
    omega = delete_rows_and_columns(omega, indices_to_remove)
    beta = beta[~np.isnan(beta)].reshape(-1, 1)
    se = se[~np.isnan(se)]
    return U, omega, beta, se


def generate_beta_se(subset_df):

    """ 
    Gathers effect sizes and standard errors from a unit of aggregation (gene).
  
    Parameters: 
    subset_df: Slice of the original dataframe that encompasses the current unit of 
        aggregation (gene).
  
    Returns: 
    beta_list: A list of effect sizes (some may be missing) from the subset.
    se_list: A list of standard errors (some may be missing) from the subset.
  
    """

    beta_list = list(subset_df["BETA"])
    se_list = list(subset_df["SE"])
    return beta_list, se_list


def calculate_all_params(
    df,
    key,
    R_phen,
    R_var_model,
    M,
    err_corr,
):

    """ 
    Calculates quantities needed for MRP (U, beta, v_beta, mu).
  
    Parameters: 
    df: Merged, filtered, and annotated dataframe containing summary statistics.
    key: gene name.
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to 
        use for analysis.
    M: Number of variants within the gene block.
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies 
        and phenotypes. Used to calculate v_beta.
  
    Returns: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K)
        dictating correlation structures, adjusted for missingness.
    beta: A S*M*K x 1 vector of effect sizes.
    v_beta: A (S*M*K x S*M*K) matrix of variances of effect sizes.
    mu: A mean of genetic effects, size of beta
        (NOTE: default is 0, can change in the code below).
  
    """

    subset_df = df[df["gene_symbol"] == key]
    sigma_m = subset_df["sigma_m_var"].tolist()
    diag_sigma_m = np.diag(np.atleast_1d(np.array(sigma_m)))
    R_var = np.diag(np.ones(M)) if R_var_model == "independent" else np.ones((M, M))
    S_var = np.dot(np.dot(diag_sigma_m, R_var), diag_sigma_m)
    beta_list, se_list = generate_beta_se(subset_df)
    beta = np.array(beta_list).reshape(-1, 1)
    se = np.array(se_list)
    omega = np.kron(err_corr, np.diag(np.ones(M)))
    R_study = np.ones((1, 1))
    U = np.kron(np.kron(R_study, R_phen), S_var)
    U, omega, beta, se = adjust_for_missingness(U, omega, beta, se, beta_list)
    diag_se = np.diag(se)
    v_beta = np.dot(np.dot(diag_se, omega), diag_se)
    mu = np.zeros(beta.shape)
    return U, beta, v_beta, mu


def output_file(bf_dfs, pop, pheno, maf_thresh, out_folder):

    """ 
    Outputs a file containing aggregation unit and Bayes Factors. 
    
    Parameters: 
    bf_dfs: List of dataframes containing Bayes Factors from each analysis.
    pop: Population name.
    pheno: Phenotype name.
    maf_thresh: Maximum MAF of variants in this run.
    out_folder: Output folder name.

    """

    outer_merge = partial(pd.merge, on="gene", how="outer")
    out_df = reduce(outer_merge, bf_dfs)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
        print("")
        print(Fore.RED + "Folder " + out_folder + " created." + Style.RESET_ALL)
        print("")
    out_file = os.path.join(
        out_folder,
        pop
        + "_"
        + pheno
        + "_"
        + str(maf_thresh)
        + ".tsv",
    )
    out_df = out_df.sort_values(by=out_df.columns[1])
    out_df.to_csv(out_file, sep="\t", index=False)
    print("")
    print(Fore.RED + "Results written to " + out_file + "." + Style.RESET_ALL)
    print("")


def get_output_file_columns(
    R_var_model,
    analysis,
):

    """
    Sets up the columns that must be enumerated in the output dataframe from MRP.

    Parameters:
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to 
        use for analysis.
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.

    Returns:
    bf_df_columns: Columns needed for the output file.
    fb: Farebrother R method (rpy2 object), or None if --p_value is not invoked.
    
    """
    if R_var_model == "independent":
        test_type = "dispersion"
    elif R_var_model == "similar":
        test_type = "burden"
    bf_df_columns = [
        "gene",
        "p_value_farebrother_" + test_type + "_" + analysis
    ]
    fb = initialize_r_objects()
    return bf_df_columns, fb


def run_mrp(
    df,
    R_phen,
    err_corr,
    R_var_model,
    analysis,
):

    """ 
    Runs MRP with the given parameters.
  
    Parameters: 
    df: Merged dataframe containing all relevant summary statistics.
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies and 
        phenotypes. Used to calculate v_beta.
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to 
        use for analysis.
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.

    Returns: 
    bf_df: Dataframe with two columns: gene and log_10 Bayes Factor. 
  
    """
    m_dict = (
        df.groupby("gene_symbol").size()
    )
    bf_df_columns, fb = get_output_file_columns(
        R_var_model,
        analysis,
    )
    data = []
    for i, (gene, value) in enumerate(m_dict.items()):
        if i % 1000 == 0:
            print("Done " + str(i) + " genes out of " + str(len(m_dict)))
        M = value
        U, beta, v_beta, mu = calculate_all_params(
            df,
            gene,
            R_phen,
            R_var_model,
            M,
            err_corr,
        )
        p_value = return_BF_pval(
            U,
            beta,
            v_beta,
            fb,
            gene,
        )
        data.append([gene] + p_value)
    bf_df = pd.DataFrame(data, columns=bf_df_columns)
    return bf_df


def print_params(
    analysis,
    R_var_model,
    maf_thresh,
):

    """ 
    Provides a text overview of each analysis in the terminal.
  
    Parameters: 
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.
        across studies.
    R_var_model: One of "independent"/"similar". Dictates correlation structure
        across variants.
    maf_thresh: Maximum MAF of variants in this run.
  
    """

    print("")
    print(Fore.YELLOW + "Analysis: " + Style.RESET_ALL + analysis)
    print(Fore.YELLOW + "R_var model: " + Style.RESET_ALL + R_var_model)
    print(Fore.YELLOW + "MAF threshold: " + Style.RESET_ALL + str(maf_thresh))
    print("")


def filter_category(df, variant_filter):

    """ 
    Filters a set of dataframes that have been read in based on functional consequence.
  
    Dependent on the variant filter that is dictated by the analysis.
  
    Parameters: 
    df: Merged dataframe containing all summary statistics.
    variant_filter: The variant filter dictated by the analysis ("ptv"/"pav"/"pcv").
  
    Returns: 
    df: Merged dataframe containing all relevant summary statistics; 
        filters out variants excluded from analysis.
  
    """

    if variant_filter == "ptv":
        df = df[df.category == "ptv"]
    elif variant_filter == "pav":
        df = df[(df.category == "ptv") | (df.category == "pav")]
    return df


def loop_through_parameters(
    df,
    maf_threshes,
    variant_filters,
    R_phen,
    R_var_models,
    err_corr,
    out_folder,
    pop,
    pheno,
):

    """ 
    Loops through parameters specified through command line (or defaults). 

    Parameters: 
    df: Merged dataframe containing all summary statistics.
    maf_threshes: List of maximum MAFs of variants in your runs.
    variant_filters: Unique list of variant filters ("ptv"/"pav"/"pcv") to use 
        for analysis.
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    R_var_models: Unique strings ("independent"/"similar") corresponding to R_var 
        matrices to use for analysis.
    err_corr: Matrix of correlation of errors across studies and phenotypes.
    out_folder: Folder where output will be placed.
  
    """

    for maf_thresh in maf_threshes:
        print(
            Fore.YELLOW
            + "Running MRP across parameters for maf_thresh "
            + str(maf_thresh)
            + "..."
            + Style.RESET_ALL
        )
        maf_df = df[(df.maf <= maf_thresh) & (df.maf > 0)]
        bf_dfs = []
        for analysis in variant_filters:
            analysis_df = filter_category(maf_df, analysis)
            for R_var_model in R_var_models:
                print_params(
                    analysis,
                    R_var_model,
                    maf_thresh,
                )
                bf_df = run_mrp(
                    analysis_df,
                    R_phen,
                    err_corr,
                    R_var_model,
                    analysis,
                )
                bf_dfs.append(bf_df)
        output_file(bf_dfs, pop, pheno, maf_thresh, out_folder)


def set_sigmas(df):

    """ 
    Assigns appropriate sigmas to appropriate variants by annotation.
  
    Filters out variants not of interest;
    Sets sigmas by functional annotation;
    Additionally adds two extra columns for two other standard choices of a uniform 
        sigma (1 and 0.05).
  
    Parameters: 
    df: Merged dataframe containing all variants across all studies and phenotypes.
  
    Returns: 
    df: Merged dataframe with an additional column (sigma_m_var) containing sigma 
        values (mapped to functional annotation via the lists inside this method).
        NOTE: One can change the sigmas associated with each type of variant by 
        adjusting the values within this method.
  
    """

    ptv = [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "start_lost",
        "stop_lost",
    ]
    pav = [
        "protein_altering_variant",
        "inframe_deletion",
        "inframe_insertion",
        "splice_region_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "missense_variant",
    ]
    proximal_coding = [
        "synonymous_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "coding_sequence_variant",
        "incomplete_terminal_codon_variant",
        "TF_binding_site_variant",
    ]
    to_filter = [
        "regulatory_region_variant",
        "intron_variant",
        "intergenic_variant",
        "downstream_gene_variant",
        "mature_miRNA_variant",
        "non_coding_transcript_exon_variant",
        "upstream_gene_variant",
        "NA",
        "NMD_transcript_variant",
    ]
    df = df[~df.most_severe_consequence.isin(to_filter)]
    sigma_m_ptv = 0.2
    sigma_m_pav = 0.05
    sigma_m_pc = 0.03
    sigma_m = dict(
        [(variant, sigma_m_ptv) for variant in ptv]
        + [(variant, sigma_m_pav) for variant in pav]
        + [(variant, sigma_m_pc) for variant in proximal_coding]
    )
    category_dict = dict(
        [(variant, "ptv") for variant in ptv]
        + [(variant, "pav") for variant in pav]
        + [(variant, "proximal_coding") for variant in proximal_coding]
    )
    sigma_m_list = list(map(sigma_m.get, df.most_severe_consequence.tolist()))
    df["sigma_m_var"] = sigma_m_list
    category_list = list(map(category_dict.get, df.most_severe_consequence.tolist()))
    df["category"] = category_list
    df = df[df.sigma_m_var.notnull()]
    return df


def merge_dfs(df, metadata_path):

    """
    Performs an outer merge on all of the files that have been read in;
    Annotates with metadata and sigma values.

    Parameters:
    df: Dataframe that contains summary statistics.
    metadata_path: Path to metadata file containing MAF, Gene symbol, etc.

    Returns:
    df: Dataframe with summary statistics and metadata.

    """

    print("")
    print(Fore.CYAN + "Merging summary statistics with metadata...")
    metadata = pd.read_csv(metadata_path, sep="\t")
    df = df.merge(metadata)
    df = set_sigmas(df)
    return df


def read_in_summary_stat(file_path):

    """
    Reads in one summary statistics file.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on 
        SE (<= 0.5).

    Parameters: 
    file_path: Path to summary statistic file.
  
    Returns: 
    df: Dataframe, ready for merge with metadata.

    """

    print(file_path)
    df = pd.read_csv(
        file_path,
        sep="\t",
        dtype={
            "#CHROM": str,
            "POS": np.int32,
            "ID": str,
            "REF": str,
            "ALT": str,
            "A1": str,
            "FIRTH?": str,
            "TEST": str,
        },
    )
    df.insert(
        loc=0,
        column="V",
        value=df["#CHROM"]
        .astype(str)
        .str.cat(df["POS"].astype(str), sep=":")
        .str.cat(df["REF"], sep=":")
        .str.cat(df["ALT"], sep=":"),
    )
    if "OR" in df.columns:
        df["BETA"] = np.log(df["OR"].astype("float64"))
    # Filter for SE as you read it in
    if "LOG(OR)_SE" in df.columns:
        df.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    df = df[df["SE"].notnull()]
    df = df[df["SE"].astype(float) <= 0.5]
    return df


def read_in_summary_stats(file_path, metadata_path, pop, pheno):

    """ 
    Reads in GBE summary statistics from the Rivas Lab file organization system on 
        Sherlock.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on 
        SE (<= 0.5).

    Contains logic for handling the case that a summary statistic file is not found.
  
    Parameters: 
    file_path: Path to summary statistic file.
    metadata_path: Path to metadata file containing MAF, Gene symbol, etc.
    pop: Population name.
    pheno: Phenotype name.
  
    Returns: 
    df: Summary statistics + metadata dataframe.

    """
    
    print("")
    print("Reading in summary statistics for:")
    print("")
    print(Fore.CYAN + "Population: " + Style.RESET_ALL + pop)
    print(Fore.CYAN + "Phenotype: " + Style.RESET_ALL + pheno)
    print("")
    try:
        df = read_in_summary_stat(file_path)
    except:
        raise IOError("File specified in --file does not exist.")
    print("")
    print(Fore.CYAN + "File found.")
    df = merge_dfs(df, metadata_path)
    return df


def print_banner():

    """ 
    Prints ASCII Art Banner + Author Info.
  
    """

    print("")
    print(Fore.RED + " __  __ ____  ____")
    print("|  \/  |  _ \|  _ \\")
    print("| |\/| | |_) | |_) |")
    print("| |  | |  _ <|  __/ ")
    print("|_|  |_|_| \_\_|  ")
    print("")
    print("Gene-based test kit" + Style.RESET_ALL)
    print("")
    print(Fore.GREEN + "Production Author:" + Style.RESET_ALL)
    print("Guhan Ram Venkataraman, B.S.H.")
    print("Ph.D. Candidate | Biomedical Informatics")
    print("")
    print(Fore.GREEN + "Contact:" + Style.RESET_ALL)
    print("Email: guhan@stanford.edu")
    print(
        "URL: https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/gene_based_test.py"
    )
    print("")
    print(Fore.GREEN + "Methods Developers:" + Style.RESET_ALL)
    print("Manuel A. Rivas, Ph.D.; Matti Pirinen, Ph.D.")
    print("Rivas Lab | Stanford University")
    print("")

def return_input_args(args):

    """ 
    Further parses the command-line input.
  
    Makes all lists unique; calculates S and K; and creates lists of appropriate 
        matrices.
  
    Parameters: 
    args: Command-line arguments that have been parsed by the parser.
  
    Returns: 
    df: Dataframe containing summary statistics and metadata.
    S: Number of populations/studies (1).
    K: Number of phenotypes (1).
    pop: Population name.
    pheno: Phenotype name.
  
    """
    pop, pheno = args.pop, args.pheno
    df = read_in_summary_stats(args.file_path, args.metadata_path, pop, pheno)
    for arg in vars(args):
        setattr(args, arg, sorted(list(set(getattr(args, arg)))))
    return (df, pop, pheno)


def range_limited_float_type(arg):

    """ 
    Type function for argparse - a float within some predefined bounds.
    
    Parameters:
    arg: Putative float.

    """

    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a valid floating point number.")
    if f <= 0 or f >= 1:
        raise argparse.ArgumentTypeError("must be > 0 and < 1.")
    return f


def initialize_parser():

    """
    Parses inputs using argparse. 
    
    Parameters: 
    valid_phenos: List of valid GBE phenotypes read in from phenotype_info.tsv.

    """

    parser = argparse.ArgumentParser(
        description="MRP takes in several variables that affect how it runs.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file",
        type=str,
        required=True,
        dest="file_path",
        help="""path to tab-separated file containing summary statistics 
       
         format of necessary columns:
        
         #CHROM   POS   REF  ALT  BETA/OR  SE
         
         """,
    )
    parser.add_argument(
        "--pop",
        type=str,
        required=True,
        dest="pop",
        help="""name of the population. used in output file naming
         """,
    )
    parser.add_argument(
        "--pheno",
        type=str,
        required=True,
        dest="pheno",
        help="""name of the phenotype. used in output file naming
         """,
    )
    parser.add_argument(
        "--metadata_path",
        type=str,
        required=True,
        dest="metadata_path",
        help="""path to tab-separated file containing:
         variants,
         gene symbols,
         consequences,
         MAFs,
         and LD independence info.
       
         format:
         
         V       gene_symbol     most_severe_consequence maf  ld_indep
         1:69081:G:C     OR4F5   5_prime_UTR_variant     0.000189471     False
        """,
    )
    parser.add_argument(
        "--R_var",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["independent"],
        dest="R_var_models",
        help="""type(s) of model across variants. 
         options: independent, similar (default: independent). can run both.
         independent is akin to the dispersion test; similar akin to burden.""",
    )
    parser.add_argument(
        "--variants",
        choices=["pcv", "pav", "ptv"],
        type=str,
        nargs="+",
        default=["ptv"],
        dest="variant_filters",
        help="""variant set(s) to consider. 
         options: proximal coding [pcv], 
                  protein-altering [pav], 
                  protein truncating [ptv] 
                  (default: ptv). can run multiple.""",
    )
    parser.add_argument(
        "--maf_thresh",
        type=range_limited_float_type,
        nargs="+",
        default=[0.01],
        dest="maf_threshes",
        help="""which MAF threshold(s) to use. must be valid floats between 0 and 1 
         (default: 0.01).""",
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        nargs=1,
        default=[],
        dest="out_folder",
        help="""folder to which output(s) will be written (default: current folder).
         if folder does not exist, it will be created.""",
    )
    return parser


if __name__ == "__main__":

    """ 
    Runs MRP analysis on GBE summary statistics with the parameters specified 
        by the command line.

    """

    import os
    parser = initialize_parser()
    args = parser.parse_args()
    print("")
    print("Valid command line arguments. Importing required packages...")
    print("")
    import pandas as pd
    from functools import partial, reduce
    pd.options.mode.chained_assignment = None
    import numpy as np
    import numpy.matlib as npm
    from numpy.linalg import LinAlgError
    from scipy.stats.stats import pearsonr
    import subprocess
    from colorama import Fore, Back, Style
    import rpy2
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.vectors import StrVector
    from rpy2.robjects.vectors import ListVector
    from rpy2.robjects.vectors import FloatVector
    import warnings
    from rpy2.rinterface import RRuntimeWarning
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    df, pop, pheno = return_input_args(args)
    out_folder = args.out_folder[0] if args.out_folder else os.getcwd()
    print_banner()
    loop_through_parameters(
        df,
        args.maf_threshes,
        args.variant_filters,
        np.ones((1, 1)),
        args.R_var_models,
        np.ones((1, 1)),
        out_folder,
        pop,
        pheno,
    )
