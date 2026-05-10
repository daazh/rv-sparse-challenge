# RV-Sparse Challenge: SpMV Implementation

This repository contains my submission for the **LFX Mentorship (RV-Sparse)**. It implements a Sparse Matrix-Vector multiplication using the **CSR (Compressed Sparse Row)** format.

## Implementation Details
- **Language**: C
- **Format**: CSR (Values, Column Indices, Row Pointers)
- **Features**: Zero-skipping logic.

## Build Instructions
```bash
gcc challenge.c -o run -lm
./run