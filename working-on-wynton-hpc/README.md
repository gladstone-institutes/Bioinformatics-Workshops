# Working on Wynton HPC

[Link to wiki](https://github.com/gladstone-institutes/Bioinformatics-Workshops/wiki/Working-on-Wynton-HPC)

## Structure

- `working-on-wynton-hpc.Rproj` : Used to manage the `renv` for the workshop slides  
  - To install the project env, open this file and run `renv::restore()` in the console   

- `./slide_materials` : Contains any images/assets needed for the slides.   
- `./workshop_materials` : Contains any materials needed by the registrants (switch to using DropBox if the size of this becomes too large)
- `./renv` : Used by `renv` to manage project library, do not modfiy these files directly  
- `./Working_on_Wynton_Part_1.Rmd` : revealjs based slides for part 1
- `./Working_on_Wynton_Part_2.Rmd` : revealjs based slides for part 2
- `./style.css` : CSS style sheet for both sets of revealjs slides
