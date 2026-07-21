---
name: Bug Report
about: Report something the Molecular AOP Analyser gets wrong
title: '[BUG] '
labels: 'bug'
assignees: 'marvinm2'
---

## Bug Description
**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

## Analysis Context
**Which part of the workflow?**
- [ ] Upload / column auto-detection
- [ ] Single analysis (`/analyze`)
- [ ] Batch wizard / comparison view
- [ ] AOP network (Cytoscape)
- [ ] Hub genes panel
- [ ] Pathway view (WikiPathways embed)
- [ ] Report generation (PDF / HTML)
- [ ] Shared results link

**Enrichment method:**
- [ ] Fisher's exact (ORA)
- [ ] GSEA
- [ ] Not applicable

**Dataset**
- Source: [uploaded file / demo dataset — which one]
- Gene ID type: [HGNC symbol / Ensembl / Entrez / other]
- Approximate number of rows:
- AOP selected: [e.g. AOP:1, NETWORK:kidney]

## Environment
- OS: [e.g., Ubuntu 20.04, macOS 12.0, Windows 10]
- Browser: [e.g., Chrome 91, Firefox 89]
- Python version: [e.g., 3.10, 3.11]
- Application version: [commit hash, image tag, or "live at molaop-analyser.vhp4safety.nl"]

**Deployment method:**
- [ ] Live instance (https://molaop-analyser.vhp4safety.nl)
- [ ] Local development (`python app.py`)
- [ ] Docker (`docker compose up`)

## Error Information
**Console errors**
If applicable, paste any JavaScript console errors:
```
Paste console errors here
```

**Server logs**
If applicable, paste relevant server log entries:
```
Paste server logs here
```

**Error messages**
Any error messages displayed to the user:
```
Paste error messages here
```

## Additional Context
**Upstream data**
- [ ] Suspected to involve KE mappings fetched from the molAOP Builder API
- [ ] Suspected to involve the AOP-Wiki RDF SPARQL endpoint
- [ ] Suspected to be a stale cache (see `cache_manager.py`)

**Additional context**
Add any other context about the problem here.

## Possible Solution
If you have a suggestion for fixing the bug, please describe it here.
