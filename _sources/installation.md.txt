# Installing Karawun

*Karawun* is python package that can be installed via *pip* or *miniconda*.

We recommend that you install *Karawun* in a python virtual environment
to avoid clashes with other packages.

There are a variety of choices of python virtual environments

## Conda based installation

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/products/individual)

1. Add package channels:
    ```bash
    conda config --append channels conda-forge --append channels anaconda --append channels SimpleITK
    ```

1. Create an environment for *Karawun*:
    ```bash
    conda create --name KarawunEnv python=3.8 karawun
    ```
    
1. Activate the environment:
    ```bash
    conda activate KarawunEnv
    ```

1. Test the installation by running the main script:
    ```bash
    importTractography -h
    ```
    This will produce the help information if the installation was successful.

Subsequent uses of karawun only require the `conda activate KarawunEnv` step.

## Pip based installation

*Karawun* is available via the [Python package index](https://pypi.org) and
can be installed using pip

1. Create and activate your choice of virtual environment.

1. Install with pip
   ```bash
   pip install karawun
   ```

1. Test the installation by running the main script:
    ```bash
    importTractography -h
    ```
    This will produce the help information if the installation was successful.

(verifying)=
## Verifying your installation

Karawun is tested before packages are distributed. However it is
recommended that you duplicate these tests locally, as follows:

1. Clone the github repository to retrieve the test data (windows users may need to install git separately)

```bash
git clone https://github.com/DevelopmentalImagingMCRI/karawun.git
```

1. Activate the environment in which Karawun is installed, and install pytest into it. eg.

```bash
conda activate KarawunEnv
pip install pytest
```

1. Change to the cloned folder

```bash
cd karawun
```

1. Run the tests

```bash
python -m pytest -s tests
```

The final line should read:

```bash
================================================== 3 passed, 2 skipped in 19.35s ===================================================
```

Skipped tests were run using custom flags to create the baseline results against
which new results are compared.