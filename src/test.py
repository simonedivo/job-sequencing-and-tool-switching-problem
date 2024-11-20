from models.SSPMF import SSPMFModel
from models.JGSMF_with_phases import solve_with_phases
from utils.test_utils import cluster_to_dataframe

cluster_directory = "../data/MTSP/Laporte/Tabela4/"

df = cluster_to_dataframe(cluster_directory)

# Get only one element for testing purposes
test_sample = df.sample(frac=1/df.shape[0], random_state=42).reset_index(drop=True)

num_jobs = test_sample['num_jobs'].values[0]
num_tools = test_sample['num_tools'].values[0]
magazine_capacity = test_sample['magazine_capacity'].values[0]
job_tools_requirements = test_sample['job_tools'].values[0]

jobs = [x for x in range(1, num_jobs + 1)]
tools = [x for x in range(1, num_tools + 1)]

SSPMF_model = SSPMFModel(jobs, tools, magazine_capacity, job_tools_requirements)
SSPMF_model.optimize()

print("Testing SSPMF model...")

SSPMF_solution = SSPMF_model.get_solution()

print("SSPMF solution")
print("Job order: ", SSPMF_solution[0], " with", SSPMF_solution[1] , "switches")

print("---------------------------------")

print("Testing JGSMF model...")

time_limit = 3600
JGSMF_solution = solve_with_phases(jobs, tools, magazine_capacity, job_tools_requirements, time_limit)

print("JGSMF solution")
print("Job order: ", JGSMF_solution[0], " with", JGSMF_solution[1] , "switches")