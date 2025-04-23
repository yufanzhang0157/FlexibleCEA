

# 🌿 FlexibleCEA


Welcome to the official repository for the research study:

**"Offering Flexible Electricity Adjustment to Power Systems Can Cut Plant Factory Energy Costs by Up to 87%"**  
by *Yufan Zhang* and *Fengqi You*, submitted to *Nature Food*.

This repository contains the code used in our study on cost reduction strategies for Controlled Environment Agriculture (CEA), focusing on optimizing electricity use in plant factories to lower energy costs.

---

## 📦 Repository Overview

This repo includes Python scripts to model, simulate, and evaluate flexible electricity use strategies in lettuce production.

### 🧩 Code Structure

| File                | Description                                                                                  |
|---------------------|----------------------------------------------------------------------------------------------|
| `economic_model.py` | Defines the optimization process for electricity and production scheduling.                  |
| `running_model.py`  | Runs the economic model over the lettuce growing cycle.                                      |
| `evaluation.py`     | Evaluates annual cost metrics, including cost and profit decomposition per kg of fresh lettuce. |

---

## ⚙️ Requirements

To run this code, please make sure the following packages are installed (versions provided for reference):

- `numpy==1.18.1`  
- `pyomo==6.6.2`

You can install them via pip:

```bash


pip install numpy==1.18.1 pyomo==6.6.2



