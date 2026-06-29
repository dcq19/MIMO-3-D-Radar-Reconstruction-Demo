# MIMO-3-D-Radar-Reconstruction-Demo

This repository provides a MATLAB demo for the first experimental dataset reported in the associated paper. The demo reconstructs a near-field 3-D image from preprocessed linear-array MIMO radar data under irregular 3-D translational motion.

For computational efficiency and easier GitHub release, the released demo uses int16 signal quantization and a modified reconstruction grid. Therefore, the numerical results may be slightly different from the full-resolution results reported in the paper.

Please download sig_data.mat from the GitHub Release page and place it in the same folder as run_mimo3d_reconstruction.m.

## Files

```text
run_mimo3d_reconstruction.m   # Main MATLAB script
params.mat                    # Parameter structure
sig_data.mat                  # Quantized pre-M2S radar signal
README.md
```

## Usage

Place all files in the same folder and run:

```matlab
run_mimo3d_reconstruction
```

The script will ask:

```text
Enable vertical motion compensation? 1 = yes, 0 = no:
```

Use `1` to enable vertical motion compensation, or `0` to disable it for comparison.

## Data

The released signal is the preprocessed MIMO echo before MIMO-to-SISO equivalent phase compensation. It has already undergone background subtraction, channel calibration, and antenna reordering.

To reduce file size, the complex signal is saved using int16 quantization. The script automatically recovers it as a single-precision complex array before reconstruction.

## Output

The script displays the reconstructed 2-D maximum-intensity projection and the measured 3-D translational motion trajectory.

## Citation

To be updated
