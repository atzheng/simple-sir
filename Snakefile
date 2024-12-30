import pandas as pd

ilinet_targets = expand(
    "results/ILINet/trajectories/{year}/{region}.csv",
    year=range(2010, 2020),
    region=range(1, 11)
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
        --samples 1   # only fit the trajectory as a whole
        """


rule sample_trajectories:
    input:
        estimates="results/ILINet/estimates/{year}/{region}.csv",
        data="data/ILINet/{year}/{region}.csv"
    output:
        "results/ILINet/trajectories/{year}/{region}.csv"
    shell:
        """
        julia --project=.. sample-trajectories.jl    \
        --estimates {input.estimates:q} \
        --data {input.data:q}           \
        --output {output:q}             \
        --seeds 20
        """


# Amazon
# -------------------------------------------------------------------------
rule preprocess_amazon:
    output:
        "data/amazon/{id}.csv"
    shell:
        """
        id={wildcards.id}            \
        input=data/amazon-raw.csv \
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
