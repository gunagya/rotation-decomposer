 
import subprocess


print('rz_argument gate_sequence')
for n in range(30):
    phase_gate_angle = f"pi/{2**n}"
    
    command_output = subprocess.check_output(f"./gridsynth {phase_gate_angle}", shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline

    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form
                               
    print(f"{n}:{phase_gate_angle} {sequence}")
    
