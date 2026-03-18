# Stochastic Eco-Evolutionary-Model (Julia)

---

## Key Features

- Trait-based population model (binned resistance levels)
- Stochastic birth, death, and mutation
- Drug response via Hill function
- Optional species interactions and pooling

## Purpose

This framework enables:
- Efficient simulation of large co-culture populations  
- Exploration of eco-evolutionary dynamics under drug pressure  
- Comparison between simulation results and experimental survival data  

## Modeling Approach

### Motivation

An exploratory, stochastic, eco-evolutionary, model in discrete-time. 
Simulates two growing, dying, mutating and interacting species under the effect of an antibiotic and regular population bottlenecks.

The objective was to construct an informed model with as few parameters as possible, while retaining the qualitative characteristics of data from evolution experiments, specifically survival times of _E. coli_ in co-culture.

Experimental populations can reach sizes of $10^7$ - $10^8$ cells, making fully explicit individual-based models (IBMs) computationally infeasible.  

To address this, I employ a **coarse-grained individual-based model**, where individuals are grouped into bins based on their resistance phenotype $R_{50}$. Each bin is updated collectively in discrete time steps, significantly reducing computational cost while retaining key evolutionary dynamics.

> **Note:**  
> $R_{50}$ is a model parameter and should not be confused with the experimental $\mathrm{IC}_{50}$.  
> - $\mathrm{IC}_{50}$: concentration where growth is reduced by half  
> - $R_{50}$: resistance level where drug-induced death is reduced by half  

---

## Model Dynamics

Let $N_{i,k}(t)$ denote the abundance of species $i$ in resistance bin $k$, with total population  
$N_i(t) = \sum_k N_{i,k}(t).$

Each time step consists of four stochastic processes applied to every bin:

### 1. Death
Cells die according to a binomial process with rate:
$d^{\mathrm{base}} + d^{\mathrm{drug}}_k(t)$

Drug-induced mortality follows a Hill function: 
$d^{\mathrm{drug}}\_{k}(t)=d^{\mathrm{drug}}\_{\mathrm{max}} \frac{c(t)^{h\_d}}{c(t)^{h\_d} + R\_{50}(k)^{h\_d}}$

---

### 2. Birth
Cell divisions are sampled from a **negative binomial distribution** to capture overdispersion at low population sizes.

Effective growth rate: $b_i^{\mathrm{base}} - b_i^{\mathrm{crowd}}(t)$

Crowding limits growth via: $b^{\mathrm{crowd}}_i(t) = \frac{K^{h_b}}{K^{h_b} + N_i(t)^{h_b}}$

- \(K = 10^8\): carrying capacity  
- Species do **not share** a carrying capacity (ensures stable coexistence across dilution cycles)

---

### 3. Mutation
- Occurs with probability $\(\mu \sim 10^{-6}\)$ per birth event  
- Mutations redistribute individuals across $\(R_{50}\)$-bins using a **log-normal distribution of fitness effects (DFE)**  
- Implemented via a precomputed mutation matrix  

---

### 4. Dilution
Every 24 hours, the population is diluted via binomial sampling with rate $\delta$.

---

## Effective Dynamics

Ignoring mutation and dilution, the expected update for each bin is: $N_{i,k}(t+\Delta t) =N_{i,k}(t)\left[1 + (b_i^{\mathrm{base}} - b_i^{\mathrm{crowd}}(t))- (d^{\mathrm{base}} + d^{\mathrm{drug}}_k(t))\right]$

---

## Mutation Establishment

After mutation, cells face stochastic extinction governed by: $p^{\mathrm{surv}}_{i,j} =1 - \frac{d^{\mathrm{eff}}_j}{b^{\mathrm{eff}}_i}$

with:
- $d^{\mathrm{eff}}_j = d^{\mathrm{base}} + d^{\mathrm{drug}}_j$
- $b^{\mathrm{eff}}_i = b_i^{\mathrm{base}}$

This creates an **effective resistance threshold**:
- Below threshold → extinction is certain  
- Above threshold → survival probability > 0  

For typical parameters:
- Threshold ≈ 2× increase in $R_{50}$  
- Maximum survival probability:
$p^{\mathrm{surv}}_{\max} \approx 0.65$

---

## Example use

_Activate the environment_
    using Pkg
    Pkg.activate("../..")
    Pkg.instantiate()

_Include all dependencies_
    include("Main Simulation.jl")

_Define the model configuration_
    co_config = ModelConfig(
        growth_fn = growth_rate,
        death_fn = death_rate,
        interaction = NoInteraction(),
        interaction_fn = growth_interaction,
        dilution = NoPooling(),
        mutation_fn = mutate_newborns,
        crowding_fn = crowd_growth,
        drug_schedule = logistic_schedule,
        record_fn = record!,
        metrics = [Every(PopSizeMetric(), 1.0)], 
        params = const_params(K = 1e7)
    )

_Run the simulation_
    results = run_simulation(co_counts0, co_config)

The results can afterwards be plotted 

    using CairoMakie
    
    fig = Figure() 
    ax = Axis(fig[1, 1],
    title = "Focal Strain Population Size over Time", 
    xlabel = "Time [hours]",
    ylabel = "Population Size",
    )
    
    CairoMakie.lines!(ax, results.PopSizeMetric.time, results.PopSizeMetric.pop1)

<img width="1200" height="900" alt="Example Plot" src="https://github.com/user-attachments/assets/0c676928-5013-43dc-ac35-a90d374aa79d" />

---

## Why this project? 

The implementation makes use of Julias dispatch feature, allowing a modular architecture and making it easy to compare different function configuration. 

In our case we modeled two different **interactions** between two strains:
1. a demographic interaction, reducing the respective population size and therefore the mutant supply.
2. an establishment interaction, suppressing the establishment of mutants.

We compared the resulting survival time distributions of the focal strain with our experimental data and were able to show that a purely demographic interaction is not sufficient to explain our results.

The modular architecture allows dynamic changes and extensions to this system, while retaining the core simulation logic.  
