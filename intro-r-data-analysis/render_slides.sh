#!/bin/bash

# Wrapper script for rendering Quarto slides to the ../docs directory
# This script renders both R Data Analysis presentations

# Set script to exit on any error
set -e

echo "Starting Quarto slide rendering..."

# Create docs directory if it doesn't exist
DOCS_DIR="../docs/Intro_to_R"
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

# Copy lesson_0 images for pre-workshop setup page
if [ -d "lesson_0/images" ]; then
    mkdir -p "$DOCS_DIR/images"
    cp -r lesson_0/images/* "$DOCS_DIR/images/"
    echo "Copied lesson_0/images to $DOCS_DIR/images"
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
render_slide "Intro_to_R_data_analysis_part_1.qmd"
render_slide "Intro_to_R_data_analysis_part_2.qmd"

# Render the pre-workshop setup page (HTML format, not revealjs)
echo "Rendering Pre-Workshop Setup page..."
quarto render "lesson_0/Pre_Workshop_Setup.qmd" --to html
if [ -f "lesson_0/Pre_Workshop_Setup.html" ]; then
    mv "lesson_0/Pre_Workshop_Setup.html" "$DOCS_DIR/"
    echo "Moved Pre_Workshop_Setup.html to $DOCS_DIR"
fi
# Clean up any supporting files
if [ -d "lesson_0/Pre_Workshop_Setup_files" ]; then
    rm -rf "lesson_0/Pre_Workshop_Setup_files"
fi

echo ""
echo "All slides rendered successfully!"
echo "Output location: $DOCS_DIR"
echo ""
echo "Generated files:"
ls -la "$DOCS_DIR"/*.html 2>/dev/null || echo "   No HTML files found"
echo ""
echo "To view the slides, open the HTML files in a web browser:"
echo "   - ${DOCS_DIR}/Intro_to_R_data_analysis_part_1.html"
echo "   - ${DOCS_DIR}/Intro_to_R_data_analysis_part_2.html"
echo "   - ${DOCS_DIR}/Pre_Workshop_Setup.html"
