 
import subprocess

PRECISION = 41 # 10^(-41)
MAX_FRACTION_POWER = 135

print('# flake8: noqa')

print('"""')
print("Gate sequences to approximate rz gates with arguments pi/2^n.")
print("Key is n so 0->rz(pi/1), 1->rz(pi/2), 2->rx(pi/4) ...")
print("Note that rz(theta) is a rotation by rx(theta/2) in our conventionprint('['")
print("")
print(f"List generated using `gridsynth` with precision -d {PRECISION}.")
print(f"It is the lowest precision for which the pi/2^{MAX_FRACTION_POWER} could be approximated.")
print('"""')
print("get_pi_over_2_to_the_n_rz_gate = [")
for n in range(MAX_FRACTION_POWER):
    phase_gate_angle = f"pi/{2**n}"
    
    command_output = subprocess.check_output(f"./gridsynth {phase_gate_angle} -d {PRECISION}", shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline

    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form
                               
    print(f"    # rz(pi/2^{n})")
    print(f"    \"{sequence}\",")

print("]")    
