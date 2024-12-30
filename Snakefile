import pandas as pd

ilinet_targets = expand(
    "results/ILINet/peaks/{type}/{year}/{region}.csv",
    year=range(2011, 2020),
    region=range(1, 11),
    type=["estimates", "estimates_true_N"]
)+ expand(
    "results/ILINet/actual_trajectories/{year}/{region}.csv",
    year=range(2011, 2020),
    region=range(1, 11),
)


amazon_products = pd.read_csv("data/amazon-products.csv").id.values
amazon_targets = expand(
    "results/Amazon/estimates/{product}.csv",
    product=amazon_products
)


rule all:
    input:
        ilinet_targets + amazon_targets


rule amazon:
    input:
        amazon_targets


rule ilinet:
    input:
        ilinet_targets



# ILINET
# -------------------------------------------------------------------------
rule download_ilinet:
    output:
        "data/ILINet.csv"
    shell:
        """
        Rscript download-ilinet.R
        """

rule preprocess_ilinet:
    input:
        "data/ILINet.csv"
    output:
        "data/ILINet/{year}/{region}.csv"
    shell:
        """
        year={wildcards.year}        \
        region={wildcards.region}    \
        output={output}              \
        input={input}                \
        j2 sql/preprocess_ilinet.sql \
        | duckdb
        """


rule fit_sir:
    input:
        "data/ILINet/{year}/{region}.csv"
    output:
        "results/ILINet/estimates/{year}/{region}.csv"
    shell:
        """
        julia --project=.. fit.jl        \
        --input {input:q}   \
        --output {output:q} \
        --model_type sir    \
        --gamma 1.68        \
        --samples 100
        """


rule fit_sir_with_true_N:
    """
    Refit the SIR model with the N estimated using all data, to demonstrate
    that estimating N is the hard part.
    """
    input:
        data="data/ILINet/{year}/{region}.csv",
        est="results/ILINet/estimates/{year}/{region}.csv"
    output:
        "results/ILINet/estimates_true_N/{year}/{region}.csv"
    shell:
        """
        julia --project=.. fit.jl \
        --input {input.data:q}    \
        --output {output:q}       \
        --model_type sir          \
        --gamma 1.68              \
        --samples 100             \
        --N $(csvcut -c N {input.est} | tail -n 1)
        """


rule get_peaks:
    """
    Compute the peak times and infection rates
    for the estimated SIR models
    """
    input:
        estimate="results/ILINet/{type}/{year}/{region}.csv",
        data="data/ILINet/{year}/{region}.csv",
    output:
        "results/ILINet/peaks/{type}/{year}/{region}.csv"
    shell:
        """
        julia --project=.. fit-flu-peaks.jl \
        --estimates {input.estimate}        \
        --data {input.data}                 \
        --output {output}                   \
        --gamma 1.68                        \
        --seeds 50
        """


rule sample_trajectories:
    """
    Sample trajectories from the estimated SIR models
    """
    input:
        estimates="results/ILINet/estimates/{year}/{region}.csv",
        data="data/ILINet/{year}/{region}.csv"
    output:
        "results/ILINet/trajectories/{year}/{region}.csv"
    shell:
        """
        julia --project=.. sample-trajectories.jl \
        --estimates {input.estimates:q}           \
        --data {input.data:q}                     \
        --output {output:q}                       \
        --seeds 20
        """


rule convert_data_to_trajectories:
    """
    Convert the ILINet data to trajectories, for visualization
    """
    input:
        "data/ILINet/{year}/{region}.csv"
    output:
        "results/ILINet/actual_trajectories/{year}/{region}.csv"
    shell:
        """
        julia --project=.. convert-data-to-trajectories.jl \
        --data {input:q} \
        --output {output:q} \
        --gamma=1.68
        """


# Amazon
# -------------------------------------------------------------------------
rule preprocess_amazon:
    output:
        "data/amazon/{id}.csv"
    shell:
        """
        id={wildcards.id}            \
        input=data/amazon-raw.csv    \
        output={output}              \
        j2 sql/preprocess_amazon.sql \
        | duckdb
        """


rule fit_bass:
    input:
        "data/amazon/{id}.csv"
    output:
        "results/Amazon/estimates/{id}.csv"
    shell:
        """
        julia --project=.. fit.jl        \
        --input {input:q}   \
        --output {output:q} \
        --model_type bass   \
        --samples 1   # only fit the trajectory as a whole
        """
