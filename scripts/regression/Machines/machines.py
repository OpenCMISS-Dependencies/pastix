import vulcain
import plafrim
import tthor
import remus
machines = {}
remote_machines = {}

for key in vulcain.machines.keys():
    machines[key] = vulcain.machines[key]
for key in vulcain.remote.keys():
    remote_machines[key] = vulcain.remote[key]

for key in remus.machines.keys():
    machines[key] = remus.machines[key]
for key in remus.remote.keys():
    remote_machines[key] = remus.remote[key]

for key in plafrim.machines.keys():
    machines[key] = plafrim.machines[key]
for key in plafrim.remote.keys():
    remote_machines[key] = plafrim.remote[key]

for key in tthor.machines.keys():
    machines[key] = tthor.machines[key]
for key in tthor.remote.keys():
    remote_machines[key] = tthor.remote[key]
