Advanced Pharmacometrics Lab - APL
Workshop materials for CaBi Data Science Camp
Date: 2025-05-30

- If you can't run the script for any reasons, open pdf versions and follow the contents for the workshop - 

1. Download current repository 
  - Click green button '<> Code' and click 'Download ZIP'
  - Extract files on your local computer

- Skip to step 5 if you have donwloaded R, Rstuio, and RTools - 

2. Download and install R
  - Go to https://www.r-project.org/ and then select the "CRAN" link under "Download"
  - Select a CRAN Mirror site - e.g. the Taiwan site is: https://cran.csie.ntu.edu.tw/
  - Select the "Download R for XXX" where XXX is your platform (e.g. Windows) 
  - If using Windows, 
      - Select the link for "install R for the first time
      - Select the link for "Previous releases" under the "Other builds" section
      - Select the link for "R 4.4.1 (June, 2024)"
  - If not using Windows, follow the links to the appropriate installs for your system

3. Download and install RStudio
  - Visit https://posit.co/download/rstudio-desktop/ to download and install the newest version of RStudio

4. Install RTools 

Windows
-------
- Download and install RTools that matches your R version (note that installing RTools will install the compiler).
- You can find RTools here: https://cran.r-project.org/bin/windows/Rtools/
- To check your R version, type `R.version` in the R console. It may be version 4.4.x or another version.

macOS
-------
- Download and install the GNU Fortran compiler to run this script on Mac.
- You can find the compiler here: https://mac.r-project.org/tools/

Linux
-------
- Ensure GCC is installed via your package manager.
    # Debian/Ubuntu
    sudo apt-get install build-essential

    # Fedora
    sudo dnf groupinstall "Development Tools"

5. Perform  test run
Once the compiler is installed
- open the CABI_DS_CAMP.Rproj file (this can be completed by double-clicking on the file or opening RStudio and opening it from there)
- within RStudio, open test.Rmd and follow the steps 
  - if errors occur (which is common), try re-running the same chunks of code
  - if errors persist, try manually performing the step where the code failed
- If you can't run the script for any reasons, open pdf files to follow the workshop

6. Open workshop materials
- After successful test run, open 2025_CABI...part1.Rmd within RStudio and follow instructions

--------------------------------

Jin Gyu Kim (kim.6182@buckeyemail.osu.edu)
