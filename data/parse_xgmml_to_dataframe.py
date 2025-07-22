import xml.etree.ElementTree as ET
import pandas as pd

# Define XGMML namespace
namespaces = {'x': 'http://www.cs.rpi.edu/XGMML'}

# Load and parse the file
xgmml_path = "wikipathways_hsa_20240410.xgmml"
tree = ET.parse(xgmml_path)
root = tree.getroot()

# Step 1: Collect a few source/target node IDs from <edge> elements
example_nodes = set()
for edge in root.findall("x:edge", namespaces):
    example_nodes.add(edge.attrib.get("source"))
    example_nodes.add(edge.attrib.get("target"))
    #if len(example_nodes) >= 5:
    #    break

# Step 2: Inspect each of those node IDs and extract all attribute key-value pairs
example_details = []

for node_id in example_nodes:
    node = next((n for n in root.findall("x:node", namespaces) if n.attrib.get("id") == node_id), None)
    if not node:
        continue

    flat_attrs = {"node_id": node_id, "label": node.attrib.get("label")}
    for att in node.findall("x:att", namespaces):
        name = att.attrib.get("name")
        value = att.attrib.get("value")
        if name:
            flat_attrs[name] = value
        for subatt in att.findall("x:att", namespaces):
            subname = subatt.attrib.get("name")
            subvalue = subatt.attrib.get("value")
            if subname:
                flat_attrs[subname] = subvalue

    example_details.append(flat_attrs)

# Step 3: Convert to DataFrame for inspection
df_node_examples = pd.DataFrame(example_details)

# Optional: display or save
# print(df_node_examples)
df_node_examples.to_csv("node_attributes.csv", index=False)

# Path to the XGMML file
xgmml_file = "wikipathways_hsa_20240410.xgmml"

# Parse the XML
tree = ET.parse(xgmml_file)
root = tree.getroot()
# Namespace handling
ns = {'ns': 'http://www.cs.rpi.edu/XGMML'}
# Extract all edge elements
edges = []
for edge in root.findall('ns:edge', ns):
    source = edge.attrib.get('source')
    target = edge.attrib.get('target')
    edge_id = edge.attrib.get('id')

    # Check if source is a WPID (starts with WP) and target is a number (gene ID)
    if source.startswith("WP") and target.isdigit():
        edges.append({
            "WPID": source,
            "gene_id": target,
            "edge_id": edge_id
        })

# Convert to DataFrame
df_edges = pd.DataFrame(edges)

print(df_edges.head())
df_edges.to_csv("edges_wpid_to_gene.csv", index=False)
