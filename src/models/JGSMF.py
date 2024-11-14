import gurobipy as gp
from gurobipy import GRB

class JGSMFModel:
    def __init(self, jobs, tools, magazine_capacity, job_tools_requirements):

        """
        jobs: list of job IDs, from 1 to N
        tools: list of tool IDs, from 1 to M
        magazine_capacity: integer, capacity of the tool magazine
        job_tools_requirements: dict, key is job ID, value is list of required tools for that job
        """

        self.jobs = jobs
        self.tools = tools
        self.magazine_capacity = magazine_capacity
        self.job_tools_requirements = job_tools_requirements # tools Ti required by job i
        self.num_nodes = list(range(0, len(jobs) + 3)) # nodes 0 and N+1 are start and end nodes, N+2 is auxiliary
        self.job_nodes = list(range(1, len(jobs) + 1)) # nodes from 1 to N
        self.N2 = len(self.num_nodes) - 1 # N+2
        self.N1 = len(self.num_nodes) - 2 # N+1

        self.model = gp.Model("JGSMF")