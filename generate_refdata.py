#!/usr/bin/env python3
"""Generate reference data XML for BoxDeformationTest.SinusoidalDeformationWorks"""

import subprocess
import sys

# Run gmx energy to extract Kinetic Energy and Potential Energy
cmd = ['./build/bin/gmx', 'energy', '-f', 'test_spc216.edr', '-o', 'energy.xvg']
input_text = "Potential\nKinetic-En.\n0\n"

try:
    result = subprocess.run(cmd, input=input_text, capture_output=True, text=True, check=True)

    # Parse energy.xvg file
    kinetic_values = []
    potential_values = []
    times = []

    with open('energy.xvg', 'r') as f:
        for line in f:
            # Skip comments and headers
            if line.startswith('#') or line.startswith('@'):
                continue

            # Parse data lines
            parts = line.split()
            if len(parts) >= 3:
                time_val = float(parts[0])
                potential_val = float(parts[1])
                kinetic_val = float(parts[2])

                times.append(time_val)
                kinetic_values.append(kinetic_val)
                potential_values.append(potential_val)

    # Calculate steps from times (dt = 0.002, nstenergy = 10)
    steps = [int(t / 0.002) for t in times]

    # Generate XML
    xml_content = '''<?xml version="1.0"?>
<?xml-stylesheet type="text/xsl" href="referencedata.xsl"?>
<ReferenceData>
  <Simulation Name="spc216">
    <Energy Name="Kinetic En.">
'''

    # Add kinetic energy values for each frame (only energy output frames)
    frame = 0
    for time, step, val in zip(times, steps, kinetic_values):
        xml_content += f'      <Real Name="Time {time:.6f} Step {step} in frame {frame}">{val:.16f}</Real>\n'
        frame += 1

    xml_content += '''    </Energy>
    <Energy Name="Potential">
'''

    # Add potential energy values
    frame = 0
    for time, step, val in zip(times, steps, potential_values):
        xml_content += f'      <Real Name="Time {time:.6f} Step {step} in frame {frame}">{val:.16f}</Real>\n'
        frame += 1

    xml_content += '''    </Energy>
  </Simulation>
</ReferenceData>
'''

    # Write to file
    output_file = 'src/programs/mdrun/tests/refdata/BoxDeformationTest_SinusoidalDeformationWorks.xml'
    with open(output_file, 'w') as f:
        f.write(xml_content)

    print(f"Generated reference data file: {output_file}")
    print(f"Found {len(times)} energy frames")
    print(f"Time range: {times[0]:.6f} - {times[-1]:.6f} ps")
    print(f"Kinetic energy range: {min(kinetic_values):.2f} - {max(kinetic_values):.2f} kJ/mol")
    print(f"Potential energy range: {min(potential_values):.2f} - {max(potential_values):.2f} kJ/mol")

except subprocess.CalledProcessError as e:
    print(f"Error running gmx energy: {e}", file=sys.stderr)
    print(f"stderr: {e.stderr}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)
