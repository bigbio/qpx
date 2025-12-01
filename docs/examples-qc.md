# Quality Control Examples

Generate comprehensive QC reports with statistics and visualizations.

## Comprehensive QC Report

Generate a complete QC report with statistics and visualizations.

```bash
#!/bin/bash

OUTPUT_DIR="output"
QC_DIR="qc_report"
mkdir -p $QC_DIR

# Generate statistics
echo "=== PSM Statistics ===" > $QC_DIR/report.txt
qpxc stats analyze psm \
    --parquet-path $OUTPUT_DIR/psm.parquet \
    >> $QC_DIR/report.txt

echo "" >> $QC_DIR/report.txt
echo "=== Project Statistics ===" >> $QC_DIR/report.txt
qpxc stats analyze project-ae \
    --absolute-path $OUTPUT_DIR/ae.parquet \
    --parquet-path $OUTPUT_DIR/psm.parquet \
    >> $QC_DIR/report.txt

# Generate QC plots
qpxc visualize plot box-intensity \
    --feature-path $OUTPUT_DIR/feature.parquet \
    --save-path $QC_DIR/intensity_boxplot.svg

qpxc visualize plot kde-intensity \
    --feature-path $OUTPUT_DIR/feature.parquet \
    --save-path $QC_DIR/intensity_kde.svg

qpxc visualize plot peptide-distribution \
    --feature-path $OUTPUT_DIR/feature.parquet \
    --save-path $QC_DIR/peptide_distribution.svg

qpxc visualize plot ibaq-distribution \
    --ibaq-path $OUTPUT_DIR/ae.parquet \
    --save-path $QC_DIR/ibaq_distribution.svg

# Create HTML report
cat > $QC_DIR/index.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>qpx QC Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; margin-top: 30px; }
        pre { background: #f4f4f4; padding: 15px; border-radius: 5px; }
        img { max-width: 100%; height: auto; margin: 20px 0; }
        .grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; }
    </style>
</head>
<body>
    <h1>qpx Quality Control Report</h1>

    <h2>Statistics</h2>
    <pre>
$(cat report.txt)
    </pre>

    <h2>Visualizations</h2>
    <div class="grid">
        <div>
            <h3>Intensity Box Plot</h3>
            <img src="intensity_boxplot.svg" alt="Intensity Box Plot">
        </div>
        <div>
            <h3>Intensity KDE</h3>
            <img src="intensity_kde.svg" alt="Intensity KDE">
        </div>
        <div>
            <h3>Peptide Distribution</h3>
            <img src="peptide_distribution.svg" alt="Peptide Distribution">
        </div>
        <div>
            <h3>iBAQ Distribution</h3>
            <img src="ibaq_distribution.svg" alt="iBAQ Distribution">
        </div>
    </div>
</body>
</html>
EOF

echo "QC report generated: $QC_DIR/index.html"
```

## Related Documentation

- **[CLI Reference](cli-reference.md)** - Complete command documentation
- **[Format Specification](format-specification.md)** - Data format details
- **[GitHub Repository](https://github.com/bigbio/qpx)** - Source code and more examples

---

**Need more examples?** Check the [`tests/examples/`](https://github.com/bigbio/qpx/tree/main/tests/examples) directory in the repository for real data files you can use for testing.

[← Back to Examples Overview](examples-overview.md)

