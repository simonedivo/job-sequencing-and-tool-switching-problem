def retrieve_format_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        first_row = lines[0].strip().split()
        num_jobs = int(first_row[0])
        num_tools = int(first_row[1])
        magazine_capacity = int(first_row[2])
        
        job_tools = {job: [] for job in range(1, num_jobs + 1)}
        
        for tool_index, line in enumerate(lines[1:], start=1):
            tool_requirements = [int(x) for x in line.strip().split()]
            for job_index, requirement in enumerate(tool_requirements, start=1):
                if requirement == 1:
                    job_tools[job_index].append(tool_index)
                    
    return num_jobs, num_tools, magazine_capacity, job_tools

def cluster_to_dataframe(cluster_directory):
    import os
    import pandas as pd

    data = []

    for filename in os.listdir(cluster_directory):
        if filename.endswith('.txt'):
            file_path = os.path.join(cluster_directory, filename)
            num_jobs, num_tools, magazine_capacity, job_tools = retrieve_format_data(file_path)
            data.append([filename, num_jobs, num_tools, magazine_capacity, job_tools])

    df = pd.DataFrame(data, columns=['filename', 'num_jobs', 'num_tools', 'magazine_capacity', 'job_tools'])
    return df