import numpy as np
import sys
import os

slurm_vars = \
{'SLURMD_NODENAME':       ('node_name',      'string'), 
 'SLURM_JOB_ID':          ('job_id',         'int'),
 'SLURM_ARRAY_JOB_ID':    ('main_job_id',    'int'),
 'SLURM_SUBMIT_DIR':      ('work_dir',       'string'),
 'SLURM_ARRAY_TASK_ID':   ('task_id',        'int'),
 'SLURM_ARRAY_TASK_STEP': ('task_step',      'int'),
 'SLURM_ARRAY_TASK_MIN':  ('task_min',       'int'),
 'SLURM_ARRAY_TASK_MAX':  ('task_max',       'int'),
 'SLURM_TASKS_PER_NODE':  ('tasks_per_node', 'int'),
 'SLURM_PROCID':          ('proc_id',        'int')}

slurm_info = dict()

for key, var in slurm_vars.items():
    name = var[0]
    cur_type = var[1]
    if key in os.environ:
        s = os.environ[key]
        if cur_type == 'string':
            slurm_info[name] = s
        elif cur_type == 'int':
            slurm_info[name] = int(s)

if 'task_id' in slurm_info:
    slurm_info['run_id'] = slurm_info['task_id'] * slurm_info['tasks_per_node'] + slurm_info['proc_id']
    slurm_info['run_num'] = (slurm_info['task_max'] + 1 - slurm_info['task_min']) // slurm_info['task_step'] * slurm_info['tasks_per_node']
else:
    slurm_info['run_id'] = 0
    slurm_info['run_num'] = 1

info = slurm_info
print(slurm_info)
