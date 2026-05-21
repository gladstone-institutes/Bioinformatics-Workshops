#!/bin/bash

# Wrapper script for rendering Quarto slides to match R Markdown output behavior
# Renders the GPU Training on CoreHPC presentation to the ../docs directory

set -e

echo "Starting Quarto slide rendering..."

DOCS_DIR="../docs/GPU_Training_on_CoreHPC"
if [ ! -d "$DOCS_DIR" ]; then
    echo "Creating output directory: $DOCS_DIR"
    mkdir -p "$DOCS_DIR"
fi

echo "Copying image assets to $DOCS_DIR..."

if [ -d "slide_materials" ]; then
    cp -r slide_materials "$DOCS_DIR/"
    echo "Copied slide_materials to $DOCS_DIR"
fi

render_slide() {
    local slide_file="$1"
    local slide_name=$(basename "$slide_file" .qmd)

    echo "Rendering $slide_file..."

    quarto render "$slide_file" --to revealjs

    if [ $? -eq 0 ]; then
        echo "Successfully rendered: ${slide_name}.html"

        if [ -f "${slide_name}.html" ]; then
            mv "${slide_name}.html" "$DOCS_DIR/"
            echo "Moved ${slide_name}.html to $DOCS_DIR"
        fi

        if [ -d "${slide_name}_files" ]; then
            rsync -av "${slide_name}_files/" "$DOCS_DIR/${slide_name}_files/"
            rm -rf "${slide_name}_files"
            echo "Synced ${slide_name}_files to $DOCS_DIR"
        fi
    else
        echo "Failed to render: $slide_file"
        exit 1
    fi
}

render_slide "GPU_Training_on_CoreHPC.qmd"

echo ""
echo "All slides rendered successfully!"
echo "Output location: $DOCS_DIR"
echo ""
echo "Generated files:"
ls -la "$DOCS_DIR"/*.html 2>/dev/null || echo "   No HTML files found"
echo ""
echo "To view the slides, open the HTML file in a web browser:"
echo "   - ${DOCS_DIR}/GPU_Training_on_CoreHPC.html"
