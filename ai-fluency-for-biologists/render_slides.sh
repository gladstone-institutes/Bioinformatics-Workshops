#!/bin/bash

# Wrapper script for rendering Quarto slides to the ../docs directory
# This script renders both AI Fluency presentations

# Set script to exit on any error
set -e

echo "Starting Quarto slide rendering..."

# Create docs directory if it doesn't exist
DOCS_DIR="../docs/AI_Fluency"
if [ ! -d "$DOCS_DIR" ]; then
    echo "Creating output directory: $DOCS_DIR"
    mkdir -p "$DOCS_DIR"
fi

# Copy image assets to docs directory
echo "Copying image assets to $DOCS_DIR..."

# Copy assets directory
if [ -d "assets" ]; then
    cp -r assets "$DOCS_DIR/"
    echo "Copied assets to $DOCS_DIR"
fi

# Function to render a single slide presentation
render_slide() {
    local slide_file="$1"
    local slide_name=$(basename "$slide_file" .qmd)

    echo "Rendering $slide_file..."

    # Render the Quarto presentation in current directory
    quarto render "$slide_file" --to revealjs

    if [ $? -eq 0 ]; then
        echo "Successfully rendered: ${slide_name}.html"

        # Move the HTML file to docs directory
        if [ -f "${slide_name}.html" ]; then
            mv "${slide_name}.html" "$DOCS_DIR/"
            echo "Moved ${slide_name}.html to $DOCS_DIR"
        fi

        # Sync the supporting files directory if it exists
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

# Render both presentations
render_slide "AI_Fluency_Part_1.qmd"
render_slide "AI_Fluency_Part_2.qmd"

echo ""
echo "All slides rendered successfully!"
echo "Output location: $DOCS_DIR"
echo ""
echo "Generated files:"
ls -la "$DOCS_DIR"/*.html 2>/dev/null || echo "   No HTML files found"
echo ""
echo "To view the slides, open the HTML files in a web browser:"
echo "   - ${DOCS_DIR}/AI_Fluency_Part_1.html"
echo "   - ${DOCS_DIR}/AI_Fluency_Part_2.html"
