# Project for the exam of Mathematical Optimisation

## Paper
The Paper on which the models are based on, is [Exploiting symmetry for the job sequencing and tool switching problem](https://www.sciencedirect.com/science/article/pii/S0377221724001632?ref=pdf_download&fr=RR-2&rr=8df4b052daf24c4a).
There is also a copy in this repository, named *paper-progetto.pdf*

## Structure
The project is structured as follows:
- `data` contains the instances needed to reproduce the experiments, in particular only instances from *data/MTSP/Laporte/Tabela4* were used
- `results` contains the two results dataframe (*new_results.csv* and *JGSMF_results.csv*) used for the scalability experiments. It also contains two scalability plots which are also visible in **scalability.ipynb**
- `src` contains:
    1. `models` which contains the two models **SSPMF** and **JGSMF**; **JGSMF** have 4 variants where each one has addeds symmetry/perfomance cuts constraints (please consider checking **scalability.ipynb** or the paper to know exactly which contraints are missing in each variant). The **JGSMF_3** is the model with all the constraints active
    2. `utils` which containts usefull functions
    3. **test.py** which reproduces a fast execution of the two models
    4. **scalability.ipynb** which containts all the scalalibity experiments

## Setup
1. Clone the repository
```sh
git clone https://github.com/simonedivo/job-sequencing-and-tool-switching-problem
```
2. Enter the repository
```sh
cd job-sequencing-and-tool-switching-problem
```
3. Install the required packages with pip or singularly by checking the *requirements.txt* file
```sh
pip install -r requirements.txt
```
4. Launch **test.py**
```sh
cd src
python test.py
```