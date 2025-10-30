# Sinusoidal Box Deformation - Custom Modification

## Overview

Custom modification to GROMACS 2025.3 adding **sinusoidal box deformation** capability for research purposes.

**Status**: ✅ Implementation complete, ⚠️ Testing required before production use
**Date**: 2025-10-30
**Base Version**: GROMACS 2025.3

---

## Motivation

Original GROMACS supports only **linear box deformation**:
```
L(t) = L(0) + rate × t
```

This modification adds **sinusoidal deformation**:
```
L(t) = L(0) + A × sin(2π × t / T)
v(t) = dL/dt = A × (2π/T) × cos(2π × t / T)
```

**Research applications:**
- Oscillatory shear simulations
- Periodic compression/expansion studies
- Material response to cyclic deformations
- Dynamic mechanical analysis (DMA) simulations

---

## Implementation Summary

### Modified Files (7 total)

| File | Purpose | Changes |
|------|---------|---------|
| `api/legacy/include/gromacs/mdtypes/md_enums.h` | Add DeformationType enum | ~15 lines |
| `api/legacy/include/gromacs/mdtypes/inputrec.h` | Add input parameters | ~5 lines |
| `src/gromacs/mdtypes/md_enums.cpp` | Enum string conversion | ~10 lines |
| `src/gromacs/gmxpreprocess/readir.cpp` | MDP parsing logic | ~75 lines |
| `src/gromacs/mdlib/boxdeformation.h` | Class interface | ~15 lines |
| `src/gromacs/mdlib/boxdeformation.cpp` | Deformation + velocity logic | ~50 lines |
| `src/gromacs/mdlib/update.cpp` | Dynamic velocity calculation | ~20 lines |

**Total**: ~190 lines added/modified

---

## Key Features

### 1. Type Selection via MDP
```mdp
deform-type = linear | sinusoidal
```

### 2. Sinusoidal Parameters
```mdp
deform-sin-amplitude = 0.5 0 0 0 0 0    ; Amplitude (nm)
deform-sin-period    = 100 0 0 0 0 0    ; Period (ps)
```
Values: `xx yy zz xy xz yz` (diagonal + shear components)

### 3. Time-Dependent Velocity Field ⚠️ CRITICAL
- **Linear deformation**: Uses constant `inputRecord.deform` (efficient)
- **Sinusoidal deformation**: Calculates `v(t) = A·ω·cos(ωt)` per timestep
- Ensures correct particle velocity adjustment for flow field

---

## Usage Examples

### Example 1: Oscillating Box Size (X-direction)
```mdp
; Box oscillates between ±0.5 nm every 100 ps
deform-type          = sinusoidal
deform-sin-amplitude = 0.5 0 0 0 0 0
deform-sin-period    = 100 0 0 0 0 0
deform-init-flow     = no
```

### Example 2: Oscillatory Shear (XY plane)
```mdp
; Shear oscillates with 1 nm amplitude, 50 ps period
deform-type          = sinusoidal
deform-sin-amplitude = 0 0 0 1.0 0 0
deform-sin-period    = 0 0 0 50 0 0
```

### Example 3: Linear Deformation (Backward Compatible)
```mdp
; Original GROMACS behavior preserved
deform-type          = linear
deform               = 0.001 0 0 0 0 0
deform-init-flow     = yes
```

---

## Technical Details

### BoxDeformation Class Changes

**New Methods:**
```cpp
void getVelocity(matrix velocity, int64_t step) const;
```
- Returns time-dependent velocity for sinusoidal deformation
- Returns constant rate for linear deformation

**apply() Method:**
```cpp
// Linear: L(t) = L(0) + rate × t
// Sinusoidal: L(t) = L(0) + A × sin(ωt)
```

### Performance Optimization

**update.cpp (line 1762-1772):**
```cpp
if (deform_ != nullptr && inputRecord.deformType == DeformationType::Sinusoidal)
{
    // Only calculate for sinusoidal (time-dependent)
    deform_->getVelocity(boxDeformationVelocity, step);
}
else
{
    // For linear, use constant inputRecord.deform (efficient)
    copy_mat(inputRecord.deform, boxDeformationVelocity);
}
```

This avoids unnecessary function calls for linear deformation.

---

## Testing Status

### ✅ Completed
- [x] Code implementation
- [x] Interface integration
- [x] Velocity field calculation
- [x] Performance optimization

### ⚠️ Required Before Production Use
- [ ] Add sinusoidal test case to `src/programs/mdrun/tests/boxdeformation.cpp`
- [ ] Validate box oscillation (amplitude, period, phase)
- [ ] Verify velocity field with analytical solution
- [ ] Energy conservation check (NVE)
- [ ] Stress tensor response
- [ ] Checkpoint/restart compatibility
- [ ] MPI parallel testing

---

## Known Limitations

1. **Single Type**: Cannot combine linear + sinusoidal simultaneously
2. **Testing**: Extensive validation needed before production
3. **Checkpointing**: May require special handling for restart
4. **Documentation**: Not in official GROMACS manual

---

## Important Warnings

⚠️ **THIS IS A CUSTOM RESEARCH MODIFICATION**

- **NOT part of official GROMACS**
- **NOT intended for GROMACS contribution**
- **Use at your own risk for research**
- **Validate results thoroughly**
- **Test with simple systems first**

---

## Future Enhancements (Optional)

If useful for research:
- [ ] Combined linear + sinusoidal
- [ ] Phase offset parameter
- [ ] Multiple frequency components
- [ ] Cosine deformation mode
- [ ] Amplitude ramp-up/down

---

## Build Instructions

```bash
cd /home/ph/programs/gromacs-2025.3
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
make -j8
make install
```

---

## Validation Checklist

Before using in production:

1. **Simple Test System**
   - Single particle in box
   - Verify box size: `L(t) = L(0) + A·sin(2πt/T)`
   - Check period and amplitude

2. **Velocity Field Test**
   - Particle at known position
   - Verify velocity correction
   - Compare with analytical flow

3. **Energy Conservation**
   - NVE ensemble
   - Monitor total energy drift
   - Should be minimal for correct implementation

4. **Stress Response**
   - Apply small amplitude oscillation
   - Measure stress tensor
   - Compare with expected viscoelastic response

---

## Citation

If using this modification in research:

**GROMACS Base:**
- Cite all standard GROMACS papers (see main README)

**This Modification:**
- Acknowledge as custom research modification
- Reference this file and modification date
- **Do not claim as official GROMACS feature**

---

## Contact & Support

This is a **research modification** with no official support.

For questions:
- Review the modified source files
- Check implementation in `boxdeformation.cpp`
- Test with simple analytical cases
- Validate results independently

---

## Version History

**2025-10-30**: Initial implementation
- Added DeformationType enum
- Implemented sinusoidal deformation logic
- Added time-dependent velocity calculation
- Performance optimization for linear case

---

**Last Updated**: 2025-10-30
**Modification Author**: Custom research implementation
**GROMACS Base Version**: 2025.3
