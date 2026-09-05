# Evolutionary Game Analysis of Coal Mine Safety Supervision in a Transition Economy: Pure Python Replication Package

This repository contains the complete 100% pure Python simulation code, parameter input datasets, and reproduction workflow for the academic paper:
> **"Phân tích trò chơi tiến hóa về giám sát an toàn mỏ than trong nền kinh tế chuyển đổi: Vai trò của chi phí thể chế, trách nhiệm chính trị và cơ chế thưởng danh tiếng"**  
> *(Evolutionary Game Analysis of Coal Mine Safety Supervision in a Transition Economy: The Roles of Institutional Costs, Political Accountability, and Reputation Incentives)*  
> Target Journal: **Safety Science / Resources Policy (Elsevier - Q1)**

---

## 1. System Requirements & Environment Setup

- **Language:** Python 3.9, 3.10, 3.11, or 3.12
- **Required Python Libraries:**
  - `numpy >= 1.22.0`
  - `scipy >= 1.8.0`
  - `matplotlib >= 3.5.0`
  - `pandas >= 1.4.0`

### Quick Installation (pip)
```bash
pip install -r requirements.txt
```

### Conda Environment Setup
```bash
conda env create -f environment.yml
conda activate coal_safety_egt
```

---

## 2. Directory Structure

```text
coal-safety-egt/
├── README.md                          # This replication guide
├── requirements.txt                   # Pip dependency specifications
├── environment.yml                    # Conda environment specifications
├── .gitignore                         # Git exclusion rules
├── main.py                            # Master Python replication script
├── data/
│   ├── administrative_benchmarks.json # Administrative coal sector data (TKV & MOLISA)
│   └── delphi_survey_data.json        # Anonymized two-round Delphi expert survey data
└── figures/                           # High-resolution vector output figures (SVG & PNG)
```

---

## 3. One-Click Full Replication Guide

To reproduce all quantitative findings, numerical verifications, Monte Carlo robustness tests, and publication-grade figures in a single step, simply run:

```bash
python main.py
```

### Modular Command-Line Execution
You can also run specific modules independently using CLI flags:
```bash
# 1. Run numerical solver verification (RK4, time-step sensitivity, and SciPy RK45 cross-check)
python main.py --verify

# 2. Compute closed-form policy thresholds (F*, H*, G*) and Jacobian eigenvalue stability
python main.py --analytical

# 3. Execute 10,000 Monte Carlo robustness iterations (reproducing Table 6)
python main.py --monte-carlo

# 4. Generate Figures 1 through 6 in both SVG and PNG format into the figures/ directory
python main.py --plots
```

---

## 4. Summary of Reproduced Results

1. **Numerical Verification & Time-Step Sensitivity:**
   - **SciPy RK45 Cross-Check:** Maximum absolute trajectory error $\max_n \|\mathbf{z}_n^{\text{RK4}} - \mathbf{z}_n^{\text{RK45}}\| < 8.5 \times 10^{-6}$; relative error $< 1.1 \times 10^{-5}$.
   - **Time-Step Sensitivity ($\Delta t = 0.05 \to 0.025$):** Full trajectory deviation $< 2.7 \times 10^{-5}$; terminal state difference at $t = 50$ is $< 1.2 \times 10^{-15}$ (change $< 0.0001\%$).
   - **Convergence Criterion:** Successive state difference $\|\mathbf{z}_{n+1} - \mathbf{z}_n\| < 10^{-7}$.

2. **Global Monte Carlo Robustness Verification ($N = 10,000$, matching Table 6):**
   - **Policy Thresholds:** $F^* = 7.56 \pm 0.82$ (95% CI $[5.97, 9.15]$), $H^* = 8.99 \pm 0.78$ (95% CI $[7.46, 10.53]$).
   - **Scenario 1 ($H = 0$):** Proportion of parameter draws yielding convergence to $E_4(1,1)$ = **0.00%** (Proving Theorem 1: Punishment alone cannot sustain the safe state).
   - **Scenario 2 ($H = 10$):** Proportion yielding convergence = **79.14%**; institutional trap $E_1(0,0)$ dismantled = **100.00%**.
   - **Scenario 3 ($H = 14$):** Proportion yielding convergence = **98.14%**; institutional trap dismantled = **100.00%**.

3. **Publication Figures (Saved in `figures/`):**
   - `fig1_phase_H0.svg` / `.png`: Phase portrait under Scenario 1 ($H = 0$, trap $E_1(0,0)$).
   - `fig2_phase_H8.svg` / `.png`: Phase portrait under Scenario 2 ($H = 10 > H^*$, safe state $E_4(1,1)$).
   - `fig3_phase_H12.svg` / `.png`: Phase portrait under Scenario 3 ($H = 14$, accelerated convergence).
   - `fig4_timeseries.svg` / `.png`: Comparative dynamic time series for $x(t)$ and $y(t)$.
   - `fig5_sensitivity_FH.svg` / `.png`: Two-factor interaction heatmap over $(F, H)$ space.
   - `fig6_sensitivity_G.svg` / `.png`: Linear marginal trade-off between subsidy $G$ and penalty $F^*$.

---

## 5. Calibrated Parameter Benchmark Table

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

---

## 6. License & Academic Citation

This simulation package is distributed under the MIT License. If you use this code or data in your academic work, please cite:
```bibtex
@article{hung_chien_2026_coal_safety,
  title={Evolutionary Game Analysis of Coal Mine Safety Supervision in a Transition Economy: The Roles of Institutional Costs, Political Accountability, and Reputation Incentives},
  author={Nguyen, Phi Hung and Trinh, Van Chien},
  journal={Safety Science / Resources Policy},
  year={2026}
}
```
