# Downloads and Plot Options

Most static plots include controls next to the figure:

- **SVG**: download a vector image generated from the current Plotly figure.
- **Options**: adjust plot width, height, and plot-specific point settings.

The options live next to each plot so sizing and export controls stay close to
the figure they affect. Filters and marker selections remain in the sidebar
because they change the data or meaning of the plot.

SVG preserves text and lines as vector elements and is suitable for reports,
slides, and downstream editing. The export uses a 1200 × 800 canvas.

## Plot size

Use **Options > Canvas** to set width and height in pixels. Larger plots are
useful when there are many facets, samples, or cell types.

## Point settings

Some plots include **Options > Points**. Depending on the plot, this can include
dot size or a jitter-dot toggle.

## Table exports

Every table includes two export controls:

- **CSV** downloads the complete underlying table as UTF-8 CSV. The browser
  grid may show only the first 2,000 rows for responsiveness, but the download
  is not truncated.
- **Excel** downloads the complete underlying table as a real `.xlsx` workbook
  with a `Data` worksheet.

Exports use the active retrieval, applied analysis settings, and current app
filters. Browser-only sorting and search filters do not change the downloaded
rows. If a table exceeds Excel's worksheet limits, use CSV instead.
