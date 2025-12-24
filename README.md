# MRes Neuroergonomics – fNIRS Data Challenge Analysis

This repository contains MATLAB scripts for the analysis of functional near-infrared spectroscopy (fNIRS) data from the MRes Neuroergonomics Data Challenge (Main Task). The analysis focuses on prefrontal cortex HbO (oxygenated hemoglobin) responses related to surgical expertise levels.

---

## Features

- Segregational analysis: comparing left vs right hemisphere activation.
- Integrational analysis: examining HbO response sensitivity to surgical expertise.
- Group comparisons across Novice, Registrar, and Consultant surgeons.
- Region of interest (ROI) analyses contrasting orbital and lateral prefrontal regions.
- Detailed visualization with box charts showcasing task-baseline differences, area under the curve (AUC), and more.
- Statistical testing including one-way ANOVA and effect size calculation (Cohen's d).

---

## Data

- The script supports loading real data from CSV files.
- Includes a simulated dataset generator for testing purposes if no real data is available.

---

## Usage

1. Modify the data loading section to point to your own fNIRS dataset in CSV format.
2. Run the full analysis script in MATLAB.
3. View statistical outputs in the MATLAB console.
4. Explore the generated figures visualizing the data by group and region.

---

## Requirements

- MATLAB R2019b or later recommended (due to `boxchart` usage).
- Statistics and Machine Learning Toolbox (for `anova1` and `multcompare`).

---

## Example Figures

- Expertise level effect on HbO Task-Baseline.
- Novice group orbital vs lateral prefrontal activation.
- Area Under Curve (AUC) by surgical expertise.

---

## Author

JUnlin Wu

---



## Contact

For questions or suggestions, please open an issue or contact [your-email@example.com].

