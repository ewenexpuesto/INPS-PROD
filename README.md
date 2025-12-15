# INPS Project — Run & Visualization Guide

Build and generate DF3 volumes

## Build the project (uses the Makefile)

```bash
make
```

## Run the main program to compute densities and export DF3 files

```bash
./main
```

This will generate files:

- `density_version0.df3` ... `density_version5.df3` are DF3 volumetric files

## Visualize with Python

```bash
python3 Plot_densities_3D.py density_version5.df3
```

## POV-Ray rendering

```bash
# Full HD with antialiasing
povray +Ipovray.pov +Odensity_df3_render_HD.png +W1920 +H1080 +A0.1 +UA

# 4K with stronger antialiasing
povray +Ipovray.pov +Odensity_df3_render_4K.png +W3840 +H2160 +A0.05 +UA
```

You might need to modify the povray. Look for "// HERE".

## Doxygen documentation

```bash
doxygen Doxyfile
```
