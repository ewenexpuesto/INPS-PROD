# INPS Project — Run & Visualization Guide

Always from the root of the project.

## Build the project

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
python3 src/Plot_densities_3D.py density_version5.df3
```

You can change the version used by changing density_version5.df3 to density_version4.df3 for instance.

## POV-Ray rendering

```bash
# Full HD with antialiasing
povray +Isrc/povray.pov +Odensity_df3_render_HD.png +W1920 +H1080 +A0.1 +UA

# 4K with stronger antialiasing
povray +Isrc/povray.pov +Odensity_df3_render_4K.png +W3840 +H2160 +A0.05 +UA
```

You might need to modify the povray. Look for "// HERE".

## Performance visualization

```bash
python3 src/Speedup_graph.py
```

## Doxygen documentation

```bash
doxygen src/Doxyfile
```
