# Fasta Viewer and Sequence Analysis Tool

This Python script provides a graphical user interface (GUI) for viewing and analyzing DNA sequences in FASTA format. The tool includes functionality to display sequence information, detect CpG islands, find homopolymers, search for motifs, and more.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)

## Requirements

- Python 3.x
- Tkinter
- pymysql

## Installation

1. Clone the repository or download the script.
2. Install the required Python packages:
    ```bash
    pip install pymysql
    ```

## Usage

Run the script using the following command:

```bash
python script_name.py
```

Replace `script_name.py` with the actual name of the script.

## Features

### Fasta Viewer GUI
- **Choose File Button**: Select a FASTA file for analysis.
- **Table Panel**: Displays sequence information in a tabular format.
- **Sequence Display**: Shows the full sequence with a scrollbar.
- **Spacer**: Displays sequence in a ruler format.
- **Motif Search**: Search for a specific motif in the sequence.
- **CpG Island Detection**: Detect CpG islands in the sequence.
- **Homopolymer Detection**: Identify homopolymers in the sequence.
- **Clear Button**: Clear all displayed data.

### Sequence Analysis
- **GC Content Calculation**: Calculates the GC content of a sequence.
- **CpG Island Detection**: Identifies CpG islands based on specified criteria.
- **Homopolymer Detection**: Finds homopolymers with a minimum length of 5.
- **Motif Search**: Locates user-defined motifs in the sequence.

### Interactive Visualization
- **Highlighting**: Highlights CpG islands in the displayed sequence.
- **Motifs**: Motifs will be lowercase when searched for in the sequence.
- **Homopolymers**: Homopolymers will be lowercase in the sequence.


