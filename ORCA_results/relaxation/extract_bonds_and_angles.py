# read the ORCA output of the relaxation
lines = open('octylsilanetriol.out').read().split('\n')

start = max(i for i, line in enumerate(lines) if 'Optimized Parameters' in line) + 3
end = next(i for i in range(start + 2, len(lines)) if lines[i].strip().startswith('---'))
table = [line for line in lines[start:end] if 'D(' not in line]

with open('bonds_and_angles.txt', 'w') as f:
    f.write('\n'.join(table) + '\n')

for number, line in enumerate(lines[start:end], start + 1):
    if 'B(' in line and 'Si' in line and 'C ' in line:
        print(f'line {number} of octylsilanetriol.out:')
        print(line)
