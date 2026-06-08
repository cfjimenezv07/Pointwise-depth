<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/Pointwise_Data_Depth_Outliers">
    <img src="Kaustlogo.png" alt="KAUST Logo" height="150">
  </a>

<h3 align="center">Pointwise data depth for univariate and multivariate functional outlier detection</h3>
</div>

## Abstract
<p align="justify">
Data depth is an efficient tool for robustly summarizing the distribution of functional data and detecting potential magnitude and shape outliers. Commonly used functional data depth notions, such as the modified band depth and extremal depth, are estimated from pointwise depth for each observed functional observation. However, these techniques require calculating one single depth value for each functional observation, which may not be sufficient to characterize the distribution of the functional data and detect potential outliers. This article presents an innovative approach to make the best use of pointwise depth. We propose using the pointwise depth distribution for magnitude outlier visualization and the correlation between pairwise depth for shape outlier detection. Furthermore, a bootstrap-based testing procedure has been introduced for the correlation to test whether there is any shape outlier. The proposed univariate methods are then extended to bivariate functional data. The performance of the proposed methods is examined and compared to conventional outlier detection techniques by intensive simulation studies. In addition, the developed methods are applied to simulated solar energy datasets from a photovoltaic system. Results revealed that the proposed method offers superior detection performance over conventional techniques. These findings will benefit engineers and practitioners in monitoring photovoltaic systems by detecting unnoticed anomalies and outliers.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Simulation Frameworks & Baseline Benchmarks (Univariate)
* **Comparisons_TVD_ED_MBD_OG.R**: Formulates the baseline univariate Monte Carlo simulation framework (500 experimental iterations, $N=1000$ curves, $M=100$ time points) across 5 distinct functional outlier contamination categories (pure dependence variance shifts, phase shift sine waves, low-amplitude high-frequency ripples, asymmetric Romo models, and locally localized central high-frequency anomalies). Quantifies and compares performance metrics (TPR and FPR) against established deep methodologies:
  * **TVD**: Total Variation Depth (Huang and Sun, 2019)
  * **ED**: Extremal Depth (Narisetty and Nair, 2016)
  * **MBD**: Modified Band Depth (López-Pintado and Romo, 2009)
  * **OG**: Outliergram (Arribas-Gil and Romo, 2014)
* **Univariate_models_MSplot.R**: Executes a parallel univariate validation study over 500 independent experimental runs ($N=1000$ curves, $M=100$ time points) to track comparative performance metrics (TPR vs. FPR) using the baseline univariate Magnitude-Shape plot framework (`fdaoutlier::msplot`).

#### 2. Simulation Frameworks & Baseline Benchmarks (Multivariate)
* **Multivariate_models_MSplot.R**: Establishes a competitive multivariate functional data benchmark loop over 500 Monte Carlo repetitions ($n=100$ multi-dimensional curves, $p=100$ evaluation steps, $d=2$ dimensions). Simulates baseline non-contaminated multivariate curves using cross-covariance bivariate Matern processes (`RandomFields::RMbiwm`) and 4 contaminated multivariate anomaly architectures to evaluate the performance of competitive Magnitude-Shape (MS) plots (`fdaoutlier::msplot`).
* **Multivariate_models_PD.R**: Implements the paper's proposed multivariate outlier discovery approach across an identical 500-iteration multivariate process. It explicitly processes multivariate points across timesteps by compounding multivariate spatial depth metrics (`fda.usc::mdepth.SD`) point-by-point to establish Pointwise Spatial Depth matrices (PSD), extracting pairwise consecutive time matrices, and measuring sample correlation structures (`cor`) to flag multivariate anomalies.

#### 3. Real-World Application & Case Studies
* **Application_PV_SI_Univ_Biv.R**: Evaluates real-world univariate and bivariate photovoltaic (PV) power systems and global solar irradiance datasets. Implements a multi-tiered outlier detection protocol by extracting magnitude outliers via functional boxplots, constructing pointwise depth (PWD) rankings, calculating sample correlations across adjacent ranks to isolate shape or dependence functional outliers, and projecting bivariate curves onto a 3D scatter plot using `scatterplot3d`.

#### 4. Replication of Figures & Visualization
* **Magnitude_PWD_Fig1.R**: Contains the foundational generation loops and graphics routines required to recreate Figure 1 of the paper. Leverages randomized Bernoulli contamination thresholds combined with symmetric Cholesky matrix factorizations over Gaussian process models to demonstrate the behavior of pointwise depths against controlled pure magnitude structural outliers.
* **PWD_PD_shape_code Fig2.R**: Computes the explicit data paths, pointwise depths (PWD), and consecutive temporal tracking steps (Pairwise Depth / PD) across 5 primary simulated models to completely recreate Figure 2 of the manuscript. Generates a multi-panel grid contrasting raw functional trajectories, boxplots of localized Sample Correlations (SC), and phase-space bivariate point distributions highlighting how contaminated trajectories shift from the central median process.

## How to Cite
If you use this code or methodology in your research, please cite our paper:
> Jiménez-Varón, C. F., Harrou, F., & Sun, Y. (2024). Pointwise data depth for univariate and multivariate functional outlier detection. *Environmetrics*, 35(5), e2851. https://doi.org/10.1002/env.2851

**BibTeX:**
```bibtex
@article{jimenezvaron2024pointwise,
  author    = {Jim{\'e}nez-Var{\'o}n, Cristian F. and Harrou, Fouzi and Sun, Ying},
  title     = {Pointwise data depth for univariate and multivariate functional outlier detection},
  journal   = {Environmetrics},
  volume    = {35},
  number    = {5},
  pages     = {e2851},
  year      = {2024},
  doi       = {10.1002/env.2851},
  url       = {[https://doi.org/10.1002/env.2851](https://doi.org/10.1002/env.2851)}
}
