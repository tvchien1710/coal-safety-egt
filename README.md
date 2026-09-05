# Evolutionary Game Analysis of Coal Mine Safety Supervision in a Transition Economy: Replication Package

Official Repository: [https://github.com/tvchien1710/coal-safety-egt](https://github.com/tvchien1710/coal-safety-egt)

This repository contains the complete pure Python simulation code, parameter input datasets, and reproduction workflow for the academic paper:
> **"Phân tích trò chơi tiến hóa về giám sát an toàn mỏ than trong nền kinh tế chuyển đổi: Vai trò của chi phí thể chế, trách nhiệm chính trị và cơ chế thưởng danh tiếng"**  
> *(Evolutionary Game Analysis of Coal Mine Safety Supervision in a Transition Economy: The Roles of Institutional Costs, Political Accountability, and Reputation Incentives)*  
> Target Journal: **Safety Science / Resources Policy (Elsevier - Q1, IF: 6.1+)**

---

## 1. System Requirements & Environment Setup

### 1.1 Python Version
- **Python >= 3.8** (Fully verified and tested on Python 3.9, 3.10, 3.11, and 3.12).

### 1.2 Library Installation
You can install dependencies using either `pip` or `conda`:

#### Option A: Pip Installation (Recommended)
```bash
pip install -r requirements.txt
```

#### Option B: Conda Environment Setup
```bash
conda env create -f environment.yml
conda activate coal_safety_egt
```

Required dependencies:
- `numpy >= 1.22.0`
- `scipy >= 1.8.0`
- `matplotlib >= 3.5.0`
- `pandas >= 1.4.0`

---

## 2. Directory Structure

```text
coal-safety-egt/
├── LICENSE                            # MIT Open Source License
├── README.md                          # Comprehensive replication documentation
├── requirements.txt                   # Pip dependency specifications
├── environment.yml                    # Conda environment specifications
├── .gitignore                         # Git exclusion rules (keeps repository clean)
├── main.py                            # Master Python CLI replication script
├── data/
│   ├── administrative_benchmarks.json # Administrative coal sector data (TKV & MOLISA)
│   └── delphi_survey_data.json        # Anonymized two-round Delphi expert survey data
└── figures/                           # Publication-grade vector (SVG) and 300 DPI (PNG) figures
    ├── fig1_phase_H0.svg / .png       # Scenario 1 phase portrait (H = 0, trap E1)
    ├── fig2_phase_H8.svg / .png       # Scenario 2 phase portrait (H = 10, safe E4)
    ├── fig3_phase_H12.svg / .png      # Scenario 3 phase portrait (H = 14, fast convergence)
    ├── fig4_timeseries.svg / .png     # Comparative dynamic time-series trajectories
    ├── fig5_sensitivity_FH.svg / .png # Two-factor interaction heatmap over (F, H) space
    └── fig6_sensitivity_G.svg / .png  # Linear marginal trade-off between G and F*
```

---

## 3. Replication Workflow & Execution Commands

### 3.1 One-Click Full Replication
To run all numerical checks, compute analytical thresholds, execute 10,000 Monte Carlo iterations, and generate all publication figures in a single command:
```bash
python main.py
```

---

### 3.2 Reproducing Table 3.1 (Numerical Error Analysis: RK4 vs RK45)
To reproduce the exact numerical error breakdown reported in **Table 3.1**:
```bash
python main.py --table3-1
```
* **Primary Metric:** Maximum absolute error $\operatorname{MaxAbsErr}_i = \max_t |z_i^{\text{RK4}}(t) - z_i^{\text{RK45}}(t)|$.
* **Supplementary Metric:** Regularized relative error $\operatorname{RelErr}_{i,t} = \frac{|z_i^{\text{RK4}}(t) - z_i^{\text{RK45}}(t)|}{|z_i^{\text{RK45}}(t)| + 10^{-4}}$, where $10^{-4}$ is a positive regularization constant preventing division by zero near boundary states ($x, y \to 0$).

Expected Table 3.1 values:
- Baseline ($z_0 = [0.5, 0.5]$):
  - $x(t)$: Max Abs Error = $3.83 \times 10^{-7}$, Max Rel Error = $5.04 \times 10^{-7}$ ($0.00005\%$)
  - $y(t)$: Max Abs Error = $1.33 \times 10^{-6}$, Max Rel Error = $1.45 \times 10^{-6}$ ($0.00015\%$)
- Full $5 \times 5$ Grid (25 Initial Conditions):
  - $x(t)$: Max Abs Error = $2.07 \times 10^{-5}$, Max Rel Error = $6.94 \times 10^{-5}$ ($< 0.007\%$)
  - $y(t)$: Max Abs Error = $2.89 \times 10^{-5}$, Max Rel Error = $6.61 \times 10^{-5}$ ($< 0.007\%$)

---

### 3.3 Reproducing Table 6 (Global Monte Carlo Robustness Verification)
To execute $10,000$ Monte Carlo iterations under $\pm 20\%$ triangular parameter uncertainty:
```bash
python main.py --monte-carlo
```
Expected Table 6 outputs:
- **Policy Thresholds (Mean $\pm$ SD [95% CI]):**
  - $F^* = 7.56 \pm 0.82 \quad [5.97, 9.15]$
  - $H^* = 8.99 \pm 0.78 \quad [7.46, 10.53]$
- **Proportion of Parameter Draws Yielding Convergence to Safe State $E_4(1, 1)$:**
  - Scenario 1 ($H = 0$): **0.00%** (Confirms Theorem 1: Fines alone cannot escape the institutional trap).
  - Scenario 2 ($H = 10 > H^*$): **79.14%** (Institutional trap dismantled in 100.00% of cases).
  - Scenario 3 ($H = 14$): **98.14%** (Accelerated convergence and high institutional resilience).

---

### 3.4 Reproducing Closed-Form Thresholds & Stability
To display analytical policy thresholds ($F^*, H^*, G^*$) and Jacobian eigenvalues ($\lambda_1, \lambda_2$) at boundary equilibria:
```bash
python main.py --analytical
```

---

### 3.5 Generating Publication Figures (Figures 1 through 6)
You can generate each figure individually or produce the entire set at once:

| Command | Generated Output | Description |
|:---|:---|:---|
| `python main.py --figure 1` | `figures/fig1_phase_H0.svg` / `.png` | **Figure 1:** Phase portrait under Scenario 1 ($H = 0$); all trajectories collapse to trap $E_1(0, 0)$. |
| `python main.py --figure 2` | `figures/fig2_phase_H8.svg` / `.png` | **Figure 2:** Phase portrait under Scenario 2 ($H = 10 > H^*$); convergence to safe state $E_4(1, 1)$. |
| `python main.py --figure 3` | `figures/fig3_phase_H12.svg` / `.png` | **Figure 3:** Phase portrait under Scenario 3 ($H = 14$); accelerated convergence to $E_4(1, 1)$. |
| `python main.py --figure 4` | `figures/fig4_timeseries.svg` / `.png` | **Figure 4:** Comparative time-series trajectories for $x(t)$ and $y(t)$ across scenarios. |
| `python main.py --figure 5` | `figures/fig5_sensitivity_FH.svg` / `.png` | **Figure 5:** Two-dimensional interaction parameter heatmap over $(F, H)$ space. |
| `python main.py --figure 6` | `figures/fig6_sensitivity_G.svg` / `.png` | **Figure 6:** Linear marginal trade-off between safety capital subsidy $G$ and minimum penalty $F^*$. |
| `python main.py --plots` | All 6 figures | Generates all figures in both vector SVG and 300 DPI PNG format. |

---

## 4. Calibrated Parameter Benchmark Table

| Parameter | Economic / Institutional Meaning | Baseline Value | Delphi Expert Range (Round 2) | Administrative Benchmark Reference |
|:---:|:---|:---:|:---:|:---|
| $p$ | Detection probability under strict inspection | **0.80** | $[0.78, 0.82]$ | Real-time multi-gas sensors, hydraulic support pressure logs |
| $q$ | Detection probability under lenient inspection | **0.20** | $[0.18, 0.22]$ | Office paperwork review, surface logbook inspections |
| $\delta$ | Ratio of lenient to strict inspection cost | **0.25** | $[0.23, 0.27]$ | Administrative overhead vs comprehensive underground audit |
| $C_r$ | Institutional cost of strict inspection | **12.0** | $[11.5, 12.6]$ | Inter-agency audit team, underground risk allowance, mobile testing |
| $C_e$ | Enterprise active safety investment cost | **3.0** | $[2.8, 3.2]$ | TKV Financial Reports: 8.5%-11.0% of production cost |
| $\Delta P$ | Excess short-term profit when bypassing safety | **4.0** | $[3.8, 4.3]$ | Short-term output surge from skipping roadway maintenance |
| $F$ | Maximum punitive fine and rectification costs | **10.0** | $[9.5, 10.6]$ | Decree 88/2020/ND-CP and Decree 12/2022/ND-CP penalties |
| $G$ | Government safety subsidy / accelerated depreciation | **1.0** | $[0.9, 1.1]$ | National Technology Innovation Fund, specialized safety tax credits |
| $L_r$ | Political accountability loss for undetected accidents | **3.0** | $[2.8, 3.2]$ | Politburo Regulation 69-QD/TW on cadre disciplinary measures |
| $H$ | Reputation incentive / promotion tournament points | **10.0** | $[9.6, 11.4]$ | Decree 90/2020/ND-CP on civil servant merit evaluation |

*Note: Probabilities $p = 0.80$ and $q = 0.20$ represent expert-informed baseline estimates elicited from 15 senior mining and inspection experts (Kendall's $W = 0.824, p < 0.001$), rather than directly observed empirical frequencies.*

---

## 5. License

This replication package and simulation software are licensed under the **MIT License**.  
See the full license terms in [LICENSE](LICENSE).

---

## 6. Citation

If you use this replication code or datasets in your research, please cite:

```bibtex
@article{hung_chien_2026_coal_safety,
  title={Evolutionary Game Analysis of Coal Mine Safety Supervision in a Transition Economy: The Roles of Institutional Costs, Political Accountability, and Reputation Incentives},
  author={Nguyen, Phi Hung and Trinh, Van Chien},
  journal={Safety Science / Resources Policy},
  year={2026},
  url={https://github.com/tvchien1710/coal-safety-egt}
}
```
