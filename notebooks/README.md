# PHYS3116 SAMI Faber–Jackson analysis

This package runs from the uploaded catalogue CSV files to the sample selection, plots, fits, comparison checks and saved numerical results. It uses numbered sections, descriptive variable names and plain-language comments.

## Open and run

1. Extract the entire ZIP. Keep the notebook, `data/` and `results/` together in the same extracted folder.
2. Open `SAMI_Faber_Jackson.ipynb` in your existing notebook environment. The included notebook already contains the outputs from a successful run.
3. Select a Python environment with NumPy, Pandas, Matplotlib and SciPy. `requirements.txt` lists those four dependencies, all used in the uploaded tutorials.
4. Restart the kernel and run the cells from the beginning. The notebook reads the provided local CSV files and regenerates the tables in `results/`.
5. `SAMI_Faber_Jackson.py` contains the same code and explanations as a Python file. Run it with the extracted package folder as the working directory. Figures appear using your Python environment's plotting backend; the notebook is the easier way to view all figures together.

Copy the whole folder into a repository subfolder such as `sami_analysis/` to keep its relative paths intact.

## What is included

| File or folder | Purpose |
|---|---|
| `SAMI_Faber_Jackson.ipynb` | Fully executed notebook with explanations and inline figures |
| `SAMI_Faber_Jackson.py` | The same analysis in a Python file |
| `data/` | Four unchanged uploaded catalogue files |
| `schemas/` | Four unchanged uploaded column-definition files |
| `results/selection_counts.csv` | Number of rows at every sequential selection stage |
| `results/observation_audit.csv` | Numerical eligibility of every original kinematic observation |
| `results/matching_audit.csv` | Whether each usable observation has the required catalogue matches |
| `results/galaxy_sample_audit.csv` | First applicable inclusion/exclusion reason for each matched galaxy |
| `results/selected_galaxies.csv` | Main sample, preserving galaxy and observation identifiers |
| `results/common_aperture_galaxies.csv` | Paired sample used for the aperture comparison |
| `results/fit_summary.csv` | Eight fits, sample sizes, coefficients, standard errors, correlations and p-values |
| `results/mass_residuals.csv` | Observed minus fitted logarithmic mass for each selected galaxy |
| `results/luminosity_exponent.csv` | r-band exponent derived from the magnitude fit |
| `results/run_log.txt` | Printed output from the checked notebook run |
| `figures/` | Seven PNG snapshots from that run |
| `METHODS.md` | Where every programming method is demonstrated in the uploaded material |
| `RESULTS.md` | Numerical findings, interpretation and remaining limits |
| `VALIDATION.md` | What was actually checked |
| `SOURCE_FILES.json` | Names, sizes and SHA-256 hashes of the supplied input copies |

## Analysis choices

The main sample requires a matching ID, positive effective-radius dispersion and error, a usable stellar mass, zero `bad_class`, and morphology code 0 or 0.5. Missing or nonpositive dispersion measurements are removed before merging. For these particular files this leaves one observation per galaxy, so no duplicate-ranking algorithm is used.

The zero-flag condition is a working selection rule. The supplied schemas identify `bad_class` as a bad/problem-object flag but do not decode the values. Confirm that definition before assessed interpretation. An explicit sensitivity fit includes the excluded nonzero-code E/E-S0 object.

The notebook uses the Sérsic effective-radius dispersion consistently for the main fit. It compares it with the 3 kpc measurement on the same galaxies and observations. It also checks strict ellipticals separately and uses spirals as a comparison population.

The primary stellar-mass fit and magnitude fit are distinct. `mstar` already contains a base-ten logarithmic mass. The logarithm is not applied twice. `m_r` is r-band absolute magnitude; no external solar-magnitude constant is introduced.

The only fitting routine is Tutorial 4's `stats.linregress`. Displayed dispersion errors are not weights in this regression. The individual standard errors do not include selection effects, uncertainty in both axes or an intrinsic-scatter model.

## Scope of the source restriction

Scientific content comes from the uploaded 2026 brief, catalogue files and schemas, tutorial material, and lecture slides. Tutorial 2: Pandas and Markdown is included in the source basis. `METHODS.md` makes this explicit. No outside paper, web source, new fitting package, or undocumented catalogue threshold was added.

Some analysis decisions and equations are applications of those taught operations, rather than verbatim tutorial exercises. The log-base conversion, residual subtraction and magnitude-slope conversion are shown algebraically. The numerical literature line is taken from the uploaded lecture slide, not from an independently retrieved paper.

The separate PNGs and `run_log.txt` are snapshots. Running the notebook updates its inline outputs and CSV files, but does not automatically refresh those snapshots or `RESULTS.md`. The analysis omits `savefig` because it was not verified in the uploaded tutorials.

This package covers programming and supporting explanations. Meeting records, attendance, past weekly commits, an individual reflection and the recorded group presentation must reflect the actual work of the group. The reference comparison is limited by the information supplied on the lecture slide.
