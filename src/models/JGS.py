import gurobipy as gp
from gurobipy import GRB

class JGSModel:
    def __init__(self, jobs, tools, magazine_capacity, job_tool_requirements):
        """
        Initialize the SSPMF model.
        
        Parameters:
        - jobs: list of job IDs.
        - tools: list of tool IDs.
        - magazine_capacity: integer, capacity of the tool magazine.
        - job_tool_requirements: dict, key is job ID, value is list of required tools for that job.
        """
        self.jobs = jobs
        self.tools = tools
        self.magazine_capacity = magazine_capacity
        self.job_tool_requirements = job_tool_requirements
        self.num_bins = len(jobs)

        self.model = gp.Model("JGS")

        self.x = None
        self.y = None
        self.v = None
        self.w = None
        self.z = None

        self.setup_variables()
        self.setup_constraints()
        self.setup_objective()

    def setup_variables(self):
        """Setup the decision variables"""

        # x[i,k]: 1 if job i is in bin k
        self.x = self.model.addVars(self.jobs, range(1, self.num_bins + 1), vtype=GRB.BINARY, name="X")

        # y[t, k]: 1 if tool t is used in bin k
        self.y = self.model.addVars(self.tools, range(1, self.num_bins + 1), vtype=gp.GRB.BINARY, name="y")
        
        # v[t, k]: 1 if tool t is loaded into bin k
        self.v = self.model.addVars(self.tools, range(1, self.num_bins + 1), vtype=gp.GRB.BINARY, name="v")
        
        # w[t, k]: 1 if tool t is removed from bin k
        self.w = self.model.addVars(self.tools, range(1, self.num_bins + 1), vtype=gp.GRB.BINARY, name="w")

        # z[k]: 1 if bin k is active
        self.z = self.model.addVars(range(1, self.num_bins + 1), vtype=gp.GRB.BINARY, name="z")

    def setup_constraints(self):
        """Define the constraints for the JGS model."""
        # Constraint: Each job must be assigned to exactly one bin
        self.model.addConstrs(
            (gp.quicksum(self.x[i, k] for k in range(1, self.num_bins + 1)) == 1 for i in self.jobs),
            name="JobAssignmentToMaxOneBin"
        )

        # Constraint: The total number of tools in each bin should not exceed magazine capacity
        self.model.addConstrs(
            (gp.quicksum(self.y[t, k] for t in self.tools) <= self.magazine_capacity for k in range(1, self.num_bins + 1)),
            name="ToolCapacityForBin"
        )

        # Constraint: For each job in a bin, all required tools for that job must be available in that bin
        for i in self.jobs:
            required_tools = self.job_tool_requirements[i]
            for t in required_tools:
                self.model.addConstrs(
                    (self.x[i, k] <= self.y[t, k] for k in range(1, self.num_bins + 1)),
                    name=f"ToolAvailability_Job{i}_Tool{t}inBin"
                )

        # Tool loading/unloading relationship constraints for tool switches
        for t in self.tools:
            for k in range(1, self.num_bins + 1):
                self.model.addConstr(
                    self.y[t, k] == self.y[t, k - 1] + self.v[t, k] - self.w[t, k],
                    name=f"ToolFlow_Tool{t}_Bin{k}"
                )

        # Simmetry Breaking Cuts

        # If a job can be done earlier, it should be done earlier
        for i in self.jobs:
            for k in range(1, self.num_bins + 1):
                self.model.addConstr(
                    self.x[i, k] <= gp.quicksum(self.x[t,k] for t in self.jobs if t < i),
                    name=f"JobBinOrdering"
                )
        
        # Bin Redundancy
        for k in range(1, self.num_bins + 1):
            self.model.addConstr(
                self.magazine_capacity * gp.quicksum(self.v[t,k] for t in self.tools) >= gp.quicksum(self.v[t,k+1] for t in self.tools),
                name=f"BinRedundancy"
            )

        # Tool Redundancy
        for k in range(1, self.num_bins + 1):
            self.model.addConstr(
                self.magazine_capacity * gp.quicksum(self.w[t,k] for t in self.tools) >= gp.quicksum(self.w[t,k+1] for t in self.tools),
                name=f"ToolRedundancy"
            )

        # Tool Redundancy 2
        for k in range(1, self.num_bins):
            for t in self.tools:
                self.model.addConstr(
                    gp.quicksum(self.x[i,k] for i in range(1, t)) >= self.w[t,k+1],
                    name=f"ToolRedundancy2"
                )

        # Job Positioning
        for i in self.jobs:
            for k in range(1, self.num_bins + 1):
                self.model.addConstr(
                    self.x[i,k] >= 1 - qp.quicksum(self.x[i,b] for b in range(1, k-1)) + gp.quicksum(self.y[t, k] for t in self.job_tool_requirements[i]) - len(self.job_tool_requirements[i]),
                    name=f"JobPositioning"
                )

        # Performance Improvement Cuts

        # Initialize z[k]
        for k in range(1, self.num_bins + 1):
            self.model.addConstr(
                self.z[k] <= gp.quicksum(self.v[t,k] for t in self.tools),
                name=f"InitializaZpart2"
            )
            for i in self.jobs:
                self.model.addConstr(
                    self.z[k] >= self.x[i,k],
                    name=f"InitializaZpart1"
                )
        
        # Alphanumeric Bin Ordering
        for k in range(1, self.num_bins):
            self.model.addConstr(
                self.z[k] >= self.z[k+1],
                name=f"AlphanumericBinOrdering"
            )
        
        # No cliques added since this model is just for testing purposes
             

    def setup_objective(self):
        """Define the objective function: minimize the number of tool switches."""
        # Objective: Minimize the total number of tool loadings (switches)
        objective = gp.quicksum(self.v[t, k] for t in self.tools for k in range(2, self.num_bins + 1))
        self.model.setObjective(objective, gp.GRB.MINIMIZE)

    def solve(self):
        """Optimize the model."""
        self.model.optimize()
        
        if self.model.status == gp.GRB.OPTIMAL:
            print("Optimal solution found.")
        elif self.model.status == gp.GRB.INFEASIBLE:
            print("No feasible solution found.")
        else:
            print("Optimization was stopped with status", self.model.status)

    def get_solution(self):
        """Retrieve the optimized job-bin assignments and tool loadings."""
        if self.model.status == gp.GRB.OPTIMAL:
            job_bin_assignments = {(i, k): self.x[i, k].x for i in self.jobs for k in range(1, self.num_bins + 1) if self.x[i, k].x > 0.5}
            tool_usage = {(t, k): self.y[t, k].x for t in self.tools for k in range(1, self.num_bins + 1) if self.y[t, k].x > 0.5}
            tool_loads = {(t, k): self.v[t, k].x for t in self.tools for k in range(2, self.num_bins + 1) if self.v[t, k].x > 0.5}
            tool_unloads = {(t, k): self.w[t, k].x for t in self.tools for k in range(2, self.num_bins + 1) if self.w[t, k].x > 0.5}
            active_bins = {k: self.z[k].x for k in range(1, self.num_bins + 1) if self.z[k].x > 0.5}
            return job_bin_assignments, tool_usage, tool_loads, tool_unloads, active_bins
        else:
            return None, None, None, None, None