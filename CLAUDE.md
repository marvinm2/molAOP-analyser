# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Application Overview

This is a Flask web application for Molecular Adverse Outcome Pathway (AOP) Key Event (KE) enrichment analysis. Users upload gene expression datasets to perform statistical enrichment testing against predefined KE gene sets, with results visualized in interactive networks using Cytoscape.js.

## Core Architecture

- **app.py**: Main Flask application with routes for file upload/preview (`/preview`) and enrichment analysis (`/analyze`)
- **helpers.py**: Contains `load_reference_sets()` function that builds KE→gene mappings from WikiPathways and AOP data
- **data/**: Contains reference datasets including AOP-KE mappings, WikiPathways data, and demo gene expression files
- **templates/**: Jinja2 templates for the web interface (index.html for upload, results.html for analysis output)
- **static/**: CSS styles and assets

## Key Data Processing Flow

1. **Reference Set Loading**: `load_reference_sets()` in helpers.py merges KE-WP mappings, WP-gene edges, and node attributes to create KE→gene dictionaries
2. **Gene Expression Processing**: User data is normalized (uppercase, duplicates combined via Fisher's method), and significance thresholds applied
3. **Enrichment Analysis**: Fisher's exact test performed for each KE using 2x2 contingency tables (significant/non-significant × in-KE/not-in-KE)
4. **Network Visualization**: Results formatted as Cytoscape.js nodes/edges for interactive AOP pathway display

## Common Development Commands

```bash
# Run locally (development)
pip install -r requirements.txt
python app.py

# Run with Docker Compose (recommended)
docker-compose up --build

# Access application
http://localhost:5000
```

### Official House Style Colors

All visualizations use the official house style color palette:

**Primary Colors:**
- Primary Dark: `#29235C` (C100 M100 Y25 K25) - Main brand color
- Primary Magenta: `#E6007E` (C0 M100 Y50 K0) - Accent color  
- Primary Blue: `#307BBF` (C80 M45 Y0 K0) - Supporting blue

**Secondary Colors:**
- Light Blue: `#009FE3` (C100 M0 Y20 K0)
- Orange: `#EB5B25` (C0 M75 Y90 K0)
- Sky Blue: `#93D5F6` (C45 M0 Y100 K0)
- Deep Magenta: `#9A1C57` (C40 M100 Y40 K10)
- Teal: `#45A6B2` (C70 M15 Y30 K0)
- Purple: `#B81178` (C30 M100 Y10 K0)
- Dark Teal: `#005A6C` (C90 M45 Y40 K30)
- Violet: `#64358C` (C75 M90 Y0 K0)

## Key Implementation Details

- Gene symbols are normalized to uppercase throughout the pipeline
- Duplicate genes (from '///' separated symbols) are combined using Fisher's method for p-values and mean for log2FC
- The application supports demo datasets (GSE90122_*.tsv files) and user uploads
- AOP networks are filtered by AOP_ID to show only relevant KEs for the selected pathway
- Enrichment results include FDR correction using Benjamini-Hochberg method
- File uploads are limited to 10MB and restricted to CSV/TSV/TXT formats

## Data Files Structure

- `aop_ke_map.csv`: Maps AOP IDs to their constituent KE IDs
- `aop_ker_edges.csv`: Defines KE relationships (edges) within each AOP
- `KE-WP.csv`: Links KE IDs to WikiPathways IDs
- `edges_wpid_to_gene.csv`: Maps WikiPathways to gene IDs
- `node_attributes.csv`: Gene ID to gene symbol mappings
- `ke_metadata.csv`: KE titles and type classifications
- always document code
- when committing to github, do not mention claude or ai
- when running the flask app, always in the background

# important-instruction-reminders
Do what has been asked; nothing more, nothing less.
NEVER create files unless they're absolutely necessary for achieving your goal.
ALWAYS prefer editing an existing file to creating a new one.
NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.
NEVER mention Claude or AI in commit messages.