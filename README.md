# funcExploreDB

This repository contains the front-end user interface for the funcExploreDB system. It connects to the locally running funcExploreDB FastAPI server to allow users to search, annotate, and explore plant stress biological entities and their relationships.

## Features

- **Live Sequence Matching**: Connects to the funcExploreDB backend (`http://localhost:8999`) to perform protein sequence searches.
- **Support for Multiple Search Engines**: Includes built-in support for searching via **ESM-C (Embed2Graph)** and **MMseqs2 (Seq2Graph)**.
- **Dynamic Relationship Extraction**: Extracts comprehensive knowledge graph relationships from the funcExploreDB database for any given protein hit or user query.
- **Adjustable Parameters**: Users can set parameters like minimum similarity threshold for ESM-C and E-value for MMseqs2 directly from the UI.
- **Sequence Hit Stacking**: Stacks redundant hits from the same plant, presenting them cleanly in a horizontally scrollable carousel.
- **Rich Visualization**: Displays evidence-supported relations and hypotheses with context (tissues, species, pathways).

## Setup & Usage

To use this front-end interface, you must first have the funcExploreDB backend server running.

1. Ensure the funcExploreDB backend server is active (typically running on `http://localhost:8999`).
2. Serve this directory using any standard HTTP server. For example:
   ```bash
   # Using Python's built-in HTTP server
   python3 -m http.server 3001
   ```
3. Open your browser and navigate to the local address (e.g., `http://localhost:3001`).

## Configuration

This demo is currently configured to point to a live backend API.

To change the API URL (e.g. for deployment), edit the `window.PSMM_API_BASE` configuration variable located at the very top of `assets/app.js`:
```javascript
// Example: Deploying to production
window.PSMM_API_BASE = 'https://api.yourdomain.com';
```

## Note on Legacy Code
Previous static data building scripts have been removed as this application now dynamically relies on the live funcExploreDB API server rather than statically built JSON bundles.
