#!/bin/bash

#SBATCH --job-name={{ job_name }}
#SBATCH --nodes={{ nodes }}
#SBATCH --ntasks-per-node={{ ntasks_per_node }}
#SBATCH --output={{ output_path }}
#SBATCH --error={{ error_path }}
#SBATCH --time={{ time }}
#SBATCH --cluster={{ cluster }}
{% if partition %}
#SBATCH --partition={{ partition }}
{% endif %}
{% if account %}
#SBATCH --account={{ account }}
{% endif %}

{% for command in commands_pre %}
command
{% endfor %}

{% for module in modules %}
module load {{ module }}
{% endfor %}

{% for key, value in env_vars.items() %}
export {{ key }}={{ value }}
{% endfor %}

{% for command in commands_run %}
command
{% endfor %}

{% for command in commands_post %}
command
{% endfor %}
