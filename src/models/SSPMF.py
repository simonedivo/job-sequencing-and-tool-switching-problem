import gurobipy as gp
from gurobipy import GRB
import math

class SSPMFModel:
    def __init__(self, jobs, tools, magazine_capacity, job_tools_requirements):

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
        
        self.model = gp.Model("SSPMF")
        self.model.setParam('OutputFlag', 0)
        self.model.setParam('TimeLimit', 3600)

        self.x = None # x[i,k]: 1 if job i is in position k (k is a node in job_nodes)
        self.y = None # y[i,j,t]: 1 if tool t is sent from node i to node j

        self.setup_variables()
        self.setup_constraints()
        self.setup_objective()

    def setup_variables(self):

        self.x = self.model.addVars(self.jobs, self.job_nodes, vtype=GRB.BINARY, name="X")
        self.y = self.model.addVars(self.num_nodes, self.num_nodes, self.tools, vtype=GRB.BINARY, name="Y")  

    def setup_constraints(self):


        for i in self.jobs:
            self.model.addConstr(gp.quicksum(self.x[i,k] for k in self.job_nodes) == 1, name=f"Job{i}InExactlyOnePosition(2b)")

        for k in self.job_nodes:
            self.model.addConstr(gp.quicksum(self.x[i,k] for i in self.jobs) == 1, name=f"EachPosition{k}HasOneJob(2c)")

        for t in self.tools:
            self.model.addConstr(self.y[0,1,t] + self.y[0,self.N2,t] == 1, name=f"Tool{t}EntersFirstNodeOrN2(2d)")

        # with range() function len(self.job_nodes)-1 is N-2 since iterate len(self.job_nodes) times but starting from 0
        for t in self.tools:
            for i in range(1, len(self.job_nodes)-1):
                self.model.addConstr(self.y[i-1,i,t] + self.y[self.N2,i,t] - self.y[i,self.N1,t] - self.y[i,i+1,t] - self.y[i,self.N2,t] == 0, name=f"FlowConservation(2e)")

        # len(self.job_nodes)-2 equal N-2 while len(self.job_nodes)-1 equals N-1, last operation should be minus like in the original paper and not plus
        for t in self.tools:
            self.model.addConstr(self.y[len(self.job_nodes)-2,len(self.job_nodes)-1,t] + self.y[self.N2,len(self.job_nodes)-1,t] - self.y[len(self.job_nodes)-1,len(self.job_nodes),t] - self.y[len(self.job_nodes)-1,self.N1,t] == 0, name=f"FlowConservationLastNode(2f)")

        for t in self.tools:
            self.model.addConstr(self.y[len(self.job_nodes)-1,len(self.job_nodes),t] - self.y[len(self.job_nodes),self.N1,t] == 0, name=f"FlowConservation(2g)")
        
        for t in self.tools:
            self.model.addConstr(gp.quicksum(self.y[i,self.N1,t] for i in self.job_nodes) == 1, name=f"FlowConservation(2h)")

        # range() of len(self.job_nodes)-1 equals N-2, there is a problem in the paper which differ from the original code paper: first sum must start at 0 not 1
        for t in self.tools:
            self.model.addConstr(gp.quicksum(self.y[i,self.N2,t] for i in range(0, len(self.job_nodes)-1)) - gp.quicksum(self.y[self.N2,i,t] for i in range(1,len(self.job_nodes))) == 0, name=f"FlowConservation(2i)")

        for i in self.job_nodes:
            required_tools = self.job_tools_requirements.get(i, [])
            for k in self.job_nodes:
                    for t in required_tools:
                        self.model.addConstr(self.x[i, k] <= self.y[k-1, k, t], name=f"UnitFlow(2j)")

        for k in self.job_nodes:
            self.model.addConstr(gp.quicksum(self.y[k-1,k,t] for t in self.tools) == self.magazine_capacity, name=f"MaxFlowCapacity(2k)")

        p = max(self.jobs, key=lambda j: len(self.job_tools_requirements.get(j, [])))
        half_sequence_length = math.ceil(len(self.job_nodes) / 2)
        self.model.addConstr(gp.quicksum(self.x[p,k] for k in range(1, half_sequence_length + 1)) == 1, name="JobMaxToolsFirstHalf(2l)")

        for t in self.tools:
            num_jobs_need_t = sum(1 for job in self.jobs if t in self.job_tools_requirements.get(job, []))       
            for k in range(1, num_jobs_need_t):
                self.model.addConstr(self.y[k, self.N1, t] == 0, name=f"Tool{t}CannotLeaveMagazineBefore{num_jobs_need_t}Jobs(2m)")

    def setup_objective(self):

        # again, with range() function the value reached is 1 unit shorter -> in range() len(self.job_nodes) is N-1 and len(self.job_nodes)-1 is N-2
        objective = gp.quicksum(self.y[i, self.N1, t] for t in self.tools for i in range(1, len(self.job_nodes))) + gp.quicksum(self.y[i, self.N2, t] for t in self.tools for i in range(1, len(self.job_nodes) - 1))
        self.model.setObjective(objective, GRB.MINIMIZE)

    def optimize(self):

        self.model.optimize()

        if self.model.status == gp.GRB.OPTIMAL:
            print("Optimal solution found.")
        elif self.model.status == gp.GRB.INFEASIBLE:
            print("No feasible solution found.")
            #self.model.computeIIS()
            #self.model.write("infeasibility_report.ilp")
        else:
            print("Optimization was stopped with status", self.model.status)

    def get_solution(self):

        if self.model.status == gp.GRB.OPTIMAL:
            job_order = sorted((k, i) for i in self.jobs for k in self.job_nodes if self.x[i, k].x > 0.5)
            job_order = [i for k, i in job_order]
            return job_order, self.count_switches()
        else:
            print("No optimal solution found.")
            return None, None
        
    def count_switches(self):
        return int(self.model.ObjVal)
        