---
phase: 06-batch-analysis
plan: GAP2
type: execute
wave: 1
depends_on: []
files_modified:
  - app.py
  - templates/batch.html
autonomous: true
requirements: [BATC-01]
gap_closure: true

must_haves:
  truths:
    - "File cards show the correct total row count (not capped at 200)"
    - "File cards show a head preview table (first 5 data rows) so users can confirm auto-detected columns against real data"
  artifacts:
    - path: "app.py"
      provides: "batch_upload() returns true total row count and head_rows list"
      contains: "sum(1 for _ in open"
    - path: "templates/batch.html"
      provides: "renderFileList() renders a head preview table per file card"
      contains: "file-card__head"
  key_links:
    - from: "app.py batch_upload()"
      to: "templates/batch.html renderFileList()"
      via: "JSON response fields row_count and head_rows"
      pattern: "head_rows.*serverFile"
---

<objective>
Fix the two UAT issues in Test 4 (File Card Metadata After Upload):

1. Row count was capped at 200 because `row_count = len(df_head)` where `df_head` uses `nrows=200`. Fix: count total file rows using a fast line count (open file, count lines, subtract 1 for header).

2. No file head preview was returned or shown. Fix: return first 5 rows of `df_head` as `head_rows` in the upload JSON response, and render a compact table in each file card so users can verify auto-detected columns against real data.

Purpose: Give users accurate row counts and visible data confirmation before committing to batch analysis settings.
Output: Fixed app.py batch_upload() + batch.html file card with head preview table.
</objective>

<execution_context>
@/home/marvin/.claude/get-shit-done/workflows/execute-plan.md
@/home/marvin/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/PROJECT.md
@.planning/ROADMAP.md
@.planning/STATE.md
@.planning/phases/06-batch-analysis/06-GAP-SUMMARY.md

<interfaces>
<!-- Key contracts for this fix. No codebase exploration needed. -->

From app.py batch_upload() — current broken behaviour (lines 1018-1039 for uploaded files, 1052-1073 for demo files):
```python
df_head = pd.read_csv(dest, sep=None, engine='python', nrows=200)
row_count = len(df_head)   # BUG: always <= 200
# ...
file_previews.append({
    'filename': safe_name,
    'row_count': row_count,
    'columns': columns,
    'suggestions': suggestions,
    # 'head_rows' missing entirely
})
```

From templates/batch.html — current JS that receives and stores upload response (lines 857-864):
```javascript
data.files.forEach(function(serverFile) {
    const entry = fileQueue.find(function(e) { return e.name === serverFile.filename; });
    if (entry) {
        entry.rowCount = serverFile.row_count;
        entry.columns = serverFile.columns;
        entry.suggestions = serverFile.suggestions;
        // entry.headRows not yet stored
    }
});
```

From templates/batch.html — current renderFileList() file card structure (lines 648-720):
The card has: icon | info(name + meta div) | removeBtn
The meta div shows: "Demo file — 1234 rows — ID: GeneSymbol, FC: log2FC, p: P.Value"
Head preview table must be appended to info div, after meta, before removeBtn.

CSS class established by 06-GAP plan:
- `.file-card__meta` — small grey text line (font-size: 0.8rem, color: #666)
New class needed:
- `.file-card__head` — compact scrollable table for head preview
</interfaces>
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix row count and add head_rows to batch_upload() response</name>
  <files>app.py</files>
  <action>
In `batch_upload()` (starting around line 985), there are two parallel blocks that each have the same bug — one for uploaded files (lines ~1018-1039) and one for demo files (lines ~1052-1073). Fix both identically:

**The fix pattern** (apply to both blocks):

Replace:
```python
df_head = pd.read_csv(dest, sep=None, engine='python', nrows=200)
row_count = len(df_head)
```

With:
```python
df_head = pd.read_csv(dest, sep=None, engine='python', nrows=200)
# Count total rows using fast line count (subtract 1 for header, handles no trailing newline)
with open(dest, 'rb') as _f:
    row_count = sum(1 for _ in _f) - 1
```

This keeps `df_head` with nrows=200 (still used for column detection and head preview), but counts real rows via byte-level line iteration — fast even for large files, no full parse.

**Add head_rows to the file_previews.append() dict** in both blocks:
```python
file_previews.append({
    'filename': safe_name,
    'row_count': row_count,
    'columns': columns,
    'suggestions': suggestions,
    'head_rows': df_head.head(5).to_dict(orient='records'),
})
```

`df_head.head(5).to_dict(orient='records')` returns a list of dicts like:
```json
[{"GeneSymbol": "ACTB", "log2FC": 0.12, "P.Value": 0.34}, ...]
```

Update the docstring for `batch_upload()` to note that `head_rows` is now included in each file dict.

Do NOT change the nrows=200 — it is still correct for column detection performance. Only fix row_count and add head_rows.
  </action>
  <verify>
Start the Flask app in the background (`python app.py &`), then upload a known file (e.g. `data/GSE90122_SR12813.tsv`) via curl or browser to `/batch/upload`. Confirm the response JSON includes `row_count` greater than 200 (this file has ~2000+ rows) and `head_rows` containing exactly 5 dicts with the expected column names.

Automated: `python -c "import ast; import subprocess; r = subprocess.run(['python', '-c', 'import app'], capture_output=True); print('Import OK' if r.returncode == 0 else r.stderr.decode())"`
  </verify>
  <done>
- `row_count` in upload response reflects actual file row count (not capped at 200)
- `head_rows` key present in each file dict in upload response with list of 5 row dicts
- Both uploaded files block and demo files block fixed identically
- App still imports without error
  </done>
</task>

<task type="auto">
  <name>Task 2: Render head preview table in batch.html file cards</name>
  <files>templates/batch.html</files>
  <action>
Two changes to batch.html:

**A. Store head_rows on entry in the upload response handler** (around line 857-864):

In the `.then(function(data) {...})` block that processes the upload response, after the existing lines that set `entry.rowCount`, `entry.columns`, `entry.suggestions`, add:
```javascript
entry.headRows = serverFile.head_rows || [];
```

**B. Add CSS for the head preview table** in the `<style>` block in `{% block head_extra %}` (place after the `.file-card__meta` rule, around line 113):

```css
/* File head preview table — shown after upload for column verification */
.file-card__head {
    margin-top: 8px;
    overflow-x: auto;
    max-height: 120px;
    overflow-y: auto;
}
.file-card__head table {
    border-collapse: collapse;
    font-size: 0.72rem;
    color: #444;
    width: max-content;
    min-width: 100%;
}
.file-card__head th,
.file-card__head td {
    border: 1px solid #ddd;
    padding: 2px 6px;
    white-space: nowrap;
    max-width: 120px;
    overflow: hidden;
    text-overflow: ellipsis;
}
.file-card__head th {
    background: #f5f5f5;
    font-weight: 600;
    position: sticky;
    top: 0;
}
```

**C. Render the head table in renderFileList()** (around line 668, inside the `fileQueue.forEach` loop):

After the block that builds and appends `meta` to `info`, and before `info.appendChild(name)` / the removeBtn section, add the head preview table. Insert it as a child of `info` after `meta`:

```javascript
// Head preview table (populated after upload)
if (entry.headRows && entry.headRows.length > 0) {
    var headWrap = document.createElement('div');
    headWrap.className = 'file-card__head';

    var tbl = document.createElement('table');
    // Header row from first entry's keys
    var thead = document.createElement('thead');
    var hrow = document.createElement('tr');
    Object.keys(entry.headRows[0]).forEach(function(col) {
        var th = document.createElement('th');
        th.textContent = col;
        hrow.appendChild(th);
    });
    thead.appendChild(hrow);
    tbl.appendChild(thead);

    // Data rows
    var tbody = document.createElement('tbody');
    entry.headRows.forEach(function(row) {
        var tr = document.createElement('tr');
        Object.values(row).forEach(function(val) {
            var td = document.createElement('td');
            td.textContent = val === null || val === undefined ? '' : String(val);
            tr.appendChild(td);
        });
        tbody.appendChild(tr);
    });
    tbl.appendChild(tbody);

    headWrap.appendChild(tbl);
    info.appendChild(headWrap);
}
```

Ensure the append order within `info` is: `name` → `meta` → `headWrap` (if present).

The head table appears only after upload completes (when `entry.headRows` is populated), so it does not affect the initial file card display before upload.
  </action>
  <verify>
Manually: Load the batch page at http://localhost:5000/batch, add a demo file (e.g. any cisplatin file), click "Next: Tag Conditions" to trigger the upload. After the upload response, the file card should show:
1. Row count reflecting the actual file size (not "200 rows")
2. A compact scrollable table showing the first 5 rows with correct column headers

Automated: `python -m py_compile templates/batch.html 2>&1 || python -c "from jinja2 import Environment, FileSystemLoader; env = Environment(loader=FileSystemLoader('templates')); env.get_template('batch.html'); print('Jinja2 parse OK')" 2>&1`
  </verify>
  <done>
- `entry.headRows` stored from upload response in fileQueue
- File cards render a `.file-card__head` table with header row + 5 data rows after upload completes
- Table uses house style colors (#f5f5f5 header background, #ddd borders)
- Table is scrollable (overflow-x: auto, max-height: 120px) to handle wide files
- Initial file card display (before upload) unchanged — no headRows table shown pre-upload
  </done>
</task>

</tasks>

<verification>
After both tasks:
1. Upload any demo file (e.g. a cisplatin dataset with thousands of rows) and confirm row count in the file card is the true file row count, not 200.
2. Confirm the file card shows a 5-row preview table with column headers matching the file's actual columns (same columns auto-populated into the column fields).
3. Confirm the existing functionality still works: column suggestions still auto-populate the id_col/fc_col/pval_col fields, tagging step still proceeds normally.
</verification>

<success_criteria>
- Row count in file cards reflects actual file row count (not capped at 200)
- File cards show a visible head preview table (5 rows, correct columns) after upload step
- No JavaScript errors in browser console during upload and render
- Existing batch upload flow (tagging, settings, analysis) remains unaffected
</success_criteria>

<output>
After completion, create `.planning/phases/06-batch-analysis/06-GAP2-SUMMARY.md` using the summary template.
</output>
